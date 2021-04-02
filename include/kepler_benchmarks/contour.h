#ifndef _KEPLER_BENCHMARKS_CONTOUR_H_
#define _KEPLER_BENCHMARKS_CONTOUR_H_

// Solve Kepler's equation via the contour integration method of Philcox et al. (2021)
// This uses techniques described in Ullisch (2020) to solve the `geometric goat problem'.

#include "math_utils.h"
#include "solver.h"

namespace kepler_benchmarks {

template <int NumGrid, typename T>
struct Contour : public Solver<T> {
  static constexpr int NumPoints = NumGrid - 2;
  T exp2R[NumPoints], exp2I[NumPoints], exp4R[NumPoints], exp4I[NumPoints], coshI[NumPoints],
      sinhI[NumPoints], ecosR[NumPoints], esinR[NumPoints];
  T radius, esinRadius, ecosRadius;

  INLINE_OR_DEVICE void precompute_for_eccentricity(const T &eccentricity) override {
    T e = this->eccentricity_ = eccentricity;

    // Define contour radius
    radius = e / 2;

    // Generate e^{ikx} sampling points and precompute real and imaginary parts
    T cf, sf, freq;
    int N_fft = (NumGrid - 1) * 2;
    for (int jj = 0; jj < NumPoints; jj++) {
      // NB: j = jj+1
      freq = 2.0 * M_PI * (jj + 1) / N_fft;
      cf = cos(freq);
      sf = sin(freq);
      exp2R[jj] = cf;
      exp2I[jj] = sf;
      exp4R[jj] = cf * cf - sf * sf;
      exp4I[jj] = 2.0 * cf * sf;
      coshI[jj] = cosh(radius * exp2I[jj]);
      sinhI[jj] = sinh(radius * exp2I[jj]);
      ecosR[jj] = e * cos(radius * exp2R[jj]);
      esinR[jj] = e * sin(radius * exp2R[jj]);
    }

    // Precompute e sin(e/2) and e cos(e/2)
    esinRadius = e * sin(radius);
    ecosRadius = e * cos(radius);
  }

  INLINE_OR_DEVICE T compute_eccentric_anomaly(const T &mean_anomaly) override {
    T this_ell = mean_anomaly;
    T e = this->eccentricity_;

    // This algorithm can't handle e == 0 so let's catch those first;
    // it also fails for M == 0; but we'll just accept that for now
    if (e < 1e-12) return mean_anomaly;

    T ft_gx2, ft_gx1, zR, zI, cosC, sinC, center;
    T fxR, fxI, ftmp, tmpcosh, tmpsinh, tmpcos, tmpsin;

    // Define contour center for each ell and precompute sin(center), cos(center)
    if (this_ell < M_PI)
      center = this_ell + e / 2;
    else
      center = this_ell - e / 2;
    sinC = sin(center);
    cosC = cos(center);

    // Accumulate Fourier coefficients
    // NB: we halve the range by symmetry, absorbing factor of 2 into ratio

    ///////////////
    // Separate out j = 0 piece, which is simpler

    // Compute z in real and imaginary parts (zI = 0 here)
    zR = center + radius;

    // Compute e*sin(zR) from precomputed quantities
    tmpsin = sinC * ecosRadius + cosC * esinRadius;  // sin(zR)

    // Compute f(z(x)) in real and imaginary parts (fxI = 0)
    fxR = zR - tmpsin - this_ell;

    // Add to array, with factor of 1/2 since an edge
    ft_gx2 = 0.5 / fxR;
    ft_gx1 = 0.5 / fxR;

    ///////////////
    // Compute for j = 1 to N_points
    // NB: j = jj+1
    for (int jj = 0; jj < NumPoints; jj++) {
      // Compute z in real and imaginary parts
      zR = center + radius * exp2R[jj];
      zI = radius * exp2I[jj];

      // Compute f(z(x)) in real and imaginary parts
      // can use precomputed cosh / sinh / cos / sin for this!
      tmpcosh = coshI[jj];                           // cosh(zI)
      tmpsinh = sinhI[jj];                           // sinh(zI)
      tmpsin = sinC * ecosR[jj] + cosC * esinR[jj];  // e sin(zR)
      tmpcos = cosC * ecosR[jj] - sinC * esinR[jj];  // e cos(zR)

      fxR = zR - tmpsin * tmpcosh - this_ell;
      fxI = zI - tmpcos * tmpsinh;

      // Compute 1/f(z) and append to array
      ftmp = fxR * fxR + fxI * fxI;
      fxR /= ftmp;
      fxI /= ftmp;

      ft_gx2 += (exp4R[jj] * fxR + exp4I[jj] * fxI);
      ft_gx1 += (exp2R[jj] * fxR + exp2I[jj] * fxI);
    }

    ///////////////
    // Separate out j = N_it piece, which is simpler

    // Compute z in real and imaginary parts (zI = 0 here)
    zR = center - radius;

    // Compute sin(zR) from precomputed quantities
    tmpsin = sinC * ecosRadius - cosC * esinRadius;  // sin(zR)

    // Compute f(z(x)) in real and imaginary parts (fxI = 0 here)
    fxR = zR - tmpsin - this_ell;

    // Add to sum, with 1/2 factor for edges
    ft_gx2 += 0.5 / fxR;
    ft_gx1 += -0.5 / fxR;

    ///////////////
    // Compute E(ell)
    return center + radius * ft_gx2 / ft_gx1;
  }

  INLINE_OR_DEVICE T compute_eccentric_anomaly(const T &mean_anomaly,
                                               const T &eccentricity) override {
    precompute_for_eccentricity(eccentricity);
    return compute_eccentric_anomaly(mean_anomaly);
  }

  INLINE_OR_DEVICE void evaluate(const T &mean_anomaly, const T &eccentricity, T *sin_ecc_anomaly,
                                 T *cos_ecc_anomaly) override {
    precompute_for_eccentricity(eccentricity);
    evaluate(mean_anomaly, sin_ecc_anomaly, cos_ecc_anomaly);
  };

  INLINE_OR_DEVICE void evaluate(const T &mean_anomaly, T *sin_ecc_anomaly,
                                 T *cos_ecc_anomaly) override {
    T E = compute_eccentric_anomaly(mean_anomaly);
    *sin_ecc_anomaly = sin(E);
    *cos_ecc_anomaly = cos(E);
  }
};

}  // namespace kepler_benchmarks

#endif