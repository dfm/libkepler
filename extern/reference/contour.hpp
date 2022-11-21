#ifndef KEPLER_EXTERN_CONTOUR_HPP
#define KEPLER_EXTERN_CONTOUR_HPP

#include <cmath>

namespace kepler {
namespace reference {

template <int NumGrid>
struct contour {
  static constexpr int NumPoints = NumGrid - 2;
  typedef double value_type;

  double eccentricity = 0.0;
  double radius, esinRadius, ecosRadius;
  double exp2R[NumPoints], exp2I[NumPoints], exp4R[NumPoints], exp4I[NumPoints], coshI[NumPoints],
      sinhI[NumPoints], ecosR[NumPoints], esinR[NumPoints];

  inline void setup(const double& eccentricity) {
    double e = this->eccentricity = eccentricity;

    // Define contour radius
    radius = e / 2;

    // Generate e^{ikx} sampling points and precompute real and imaginary parts
    double cf, sf, freq;
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

  inline double solve(const double& mean_anomaly) const {
    double this_ell = mean_anomaly;
    double e = this->eccentricity;

    // This algorithm can't handle e == 0 or M == 0 so let's catch those first
    if (e < 1e-12 || this_ell < 1e-12) return mean_anomaly;

    double ft_gx2, ft_gx1, zR, zI, cosC, sinC, center;
    double fxR, fxI, ftmp, tmpcosh, tmpsinh, tmpcos, tmpsin;

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
};

}  // namespace reference
}  // namespace kepler

#endif
