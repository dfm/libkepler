#ifndef _KEPLER_BENCHMARKS_BRANDT_H_
#define _KEPLER_BENCHMARKS_BRANDT_H_

// C code to compute the eccentric anomaly, its sine and cosine, from
// an input mean anomaly and eccentricity.  Timothy D. Brandt wrote
// the code; the algorithm is based on Raposo-Pulido & Pelaez, 2017,
// MNRAS, 467, 1702.  Computational cost is equivalent to around 3
// trig calls (sine, cosine) in tests as of August 2020.  This can be
// further reduced if using many mean anomalies at fixed eccentricity;
// the code would need some modest refactoring in that case.  Accuracy
// should be within a factor of a few of machine epsilon in E-ecc*sinE
// and in cosE up to at least ecc=0.999999.

#include <cstdlib>

#include "math_utils.h"
#include "solver.h"

namespace kepler_benchmarks {
namespace brandt {

// Evaluate sine with a series expansion.  We can guarantee that the
// argument will be <=pi/4, and this reaches double precision (within
// a few machine epsilon) at a significantly lower cost than the
// function call to sine that obeys the IEEE standard.
template <typename Scalar>
INLINE_OR_DEVICE Scalar shortsin(const Scalar &x) {
  const Scalar if3 = 1. / 6;
  const Scalar if5 = 1. / (6. * 20);
  const Scalar if7 = 1. / (6. * 20 * 42);
  const Scalar if9 = 1. / (6. * 20 * 42 * 72);
  const Scalar if11 = 1. / (6. * 20 * 42 * 72 * 110);
  const Scalar if13 = 1. / (6. * 20 * 42 * 72 * 110 * 156);
  const Scalar if15 = 1. / (6. * 20 * 42 * 72 * 110 * 156 * 210);

  Scalar x2 = x * x;
  return x *
         (1 - x2 * (if3 -
                    x2 * (if5 - x2 * (if7 - x2 * (if9 - x2 * (if11 - x2 * (if13 - x2 * if15)))))));
}

// Use the second-order series expanion in Raposo-Pulido & Pelaez
// (2017) in the singular corner (eccentricity close to 1, mean
// anomaly close to zero).
template <typename Scalar>
INLINE_OR_DEVICE Scalar EAstart(const Scalar &M, const Scalar &ecc) {
  const Scalar ome = 1. - ecc;
  const Scalar sqrt_ome = sqrt(ome);

  const Scalar chi = M / (sqrt_ome * ome);
  const Scalar Lam = sqrt(8 + 9 * chi * chi);
  const Scalar S = cbrt(Lam + 3 * chi);
  const Scalar sigma = 6 * chi / (2 + S * S + 4. / (S * S));
  const Scalar s2 = sigma * sigma;
  const Scalar s4 = s2 * s2;

  const Scalar denom = 1.0 / (s2 + 2);
  const Scalar E =
      sigma * (1 + s2 * ome * denom *
                       ((s2 + 20) / 60. +
                        s2 * ome * denom * denom * (s2 * s4 + 25 * s4 + 340 * s2 + 840) / 1400));

  return E * sqrt_ome;
}

// Calculate the eccentric anomaly, its sine and cosine, using a
// variant of the algorithm suggested in Raposo-Pulido & Pelaez (2017)
// and used in Brandt et al. (2020).  Use the series expansion above
// to generate an initial guess in the singular corner and use a
// fifth-order polynomial to get the initial guess otherwise.  Use
// series and square root calls to evaluate sine and cosine, and
// update their values using series.  Accurate to better than 1e-15 in
// E-ecc*sin(E)-M at all mean anomalies and at eccentricies up to
// 0.999999.
template <typename Scalar>
INLINE_OR_DEVICE Scalar calcEA(const Scalar &M, const Scalar &ecc, Scalar *sinE, Scalar *cosE) {
  if (ecc < 1e-12 || M < 1e-12) {
    *sinE = 0;
    *cosE = 1;
    return M;
  }

  const Scalar one_sixth = 1. / 6;

  const Scalar pi = M_PI;
  const Scalar pi_d_12 = M_PI / 12;
  const Scalar pi_d_6 = M_PI / 6;
  const Scalar pi_d_4 = M_PI / 4;
  const Scalar pi_d_3 = M_PI / 3;
  const Scalar fivepi_d_12 = M_PI * 5. / 12;
  const Scalar pi_d_2 = M_PI / 2;
  const Scalar sevenpi_d_12 = M_PI * 7. / 12;
  const Scalar twopi_d_3 = M_PI * 2. / 3;
  const Scalar threepi_d_4 = M_PI * 3. / 4;
  const Scalar fivepi_d_6 = M_PI * 5. / 6;
  const Scalar elevenpi_d_12 = M_PI * 11. / 12;
  const Scalar twopi = M_PI * 2;

  const Scalar g2s_e = 0.2588190451025207623489 * ecc;
  const Scalar g3s_e = 0.5 * ecc;
  const Scalar g4s_e = 0.7071067811865475244008 * ecc;
  const Scalar g5s_e = 0.8660254037844386467637 * ecc;
  const Scalar g6s_e = 0.9659258262890682867497 * ecc;

  Scalar bounds[13];
  Scalar EA_tab[9];

  int k;
  Scalar MA, EA, sE, cE, x, y;
  Scalar B0, B1, B2, dx, idx;
  int MAsign = 1;
  Scalar one_over_ecc = 1e17;
  if (ecc > 1e-17) one_over_ecc = 1. / ecc;

  MA = M;
  // MA = MAmod(M);
  if (MA > pi) {
    MAsign = -1;
    MA = twopi - MA;
  }

  // Series expansion
  if (2 * MA + 1 - ecc < 0.2) {
    EA = EAstart(MA, ecc);
  } else {
    // Polynomial boundaries given in Raposo-Pulido & Pelaez
    bounds[0] = 0;
    bounds[1] = pi_d_12 - g2s_e;
    bounds[2] = pi_d_6 - g3s_e;
    bounds[3] = pi_d_4 - g4s_e;
    bounds[4] = pi_d_3 - g5s_e;
    bounds[5] = fivepi_d_12 - g6s_e;
    bounds[6] = pi_d_2 - ecc;
    bounds[7] = sevenpi_d_12 - g6s_e;
    bounds[8] = twopi_d_3 - g5s_e;
    bounds[9] = threepi_d_4 - g4s_e;
    bounds[10] = fivepi_d_6 - g3s_e;
    bounds[11] = elevenpi_d_12 - g2s_e;
    bounds[12] = pi;

    // Which interval?
    for (k = 11; k > 0; k--) {
      if (MA > bounds[k]) break;
    }
    // if (k < 0) k = 0;

    // Values at the two endpoints.
    EA_tab[0] = k * pi_d_12;
    EA_tab[6] = (k + 1) * pi_d_12;

    // First two derivatives at the endpoints. Left endpoint first.
    int sign = (k >= 6) ? 1 : -1;

    x = 1 / (1 - ((6 - k) * pi_d_12 + sign * bounds[abs(6 - k)]));
    y = -0.5 * (k * pi_d_12 - bounds[k]);
    EA_tab[1] = x;
    EA_tab[2] = y * x * x * x;

    x = 1 / (1 - ((5 - k) * pi_d_12 + sign * bounds[abs(5 - k)]));
    y = -0.5 * ((k + 1) * pi_d_12 - bounds[k + 1]);
    EA_tab[7] = x;
    EA_tab[8] = y * x * x * x;

    // Solve a matrix equation to get the rest of the coefficients.
    idx = 1 / (bounds[k + 1] - bounds[k]);

    B0 = idx * (-EA_tab[2] - idx * (EA_tab[1] - idx * pi_d_12));
    B1 = idx * (-2 * EA_tab[2] - idx * (EA_tab[1] - EA_tab[7]));
    B2 = idx * (EA_tab[8] - EA_tab[2]);

    EA_tab[3] = B2 - 4 * B1 + 10 * B0;
    EA_tab[4] = (-2 * B2 + 7 * B1 - 15 * B0) * idx;
    EA_tab[5] = (B2 - 3 * B1 + 6 * B0) * idx * idx;

    // Now use the coefficients of this polynomial to get the initial guess.
    dx = MA - bounds[k];
    EA =
        EA_tab[0] +
        dx * (EA_tab[1] + dx * (EA_tab[2] + dx * (EA_tab[3] + dx * (EA_tab[4] + dx * EA_tab[5]))));
  }

  // Sine and cosine initial guesses using series
  if (EA < pi_d_4) {
    sE = shortsin(EA);
    cE = sqrt(1 - sE * sE);
  } else if (EA > threepi_d_4) {
    sE = shortsin(pi - EA);
    cE = -sqrt(1 - sE * sE);
  } else {
    cE = shortsin(pi_d_2 - EA);
    sE = sqrt(1 - cE * cE);
  }

  Scalar num, denom, dEA;

  // Halley's method to update E
  num = (MA - EA) * one_over_ecc + sE;
  denom = one_over_ecc - cE;
  dEA = num * denom / (denom * denom + 0.5 * sE * num);

  // Use series to update sin and cos
  Scalar E = 0;
  if (ecc < 0.78 || MA > 0.4) {
    E = MAsign * (EA + dEA);
    *sinE = MAsign * (sE * (1 - 0.5 * dEA * dEA) + dEA * cE);
    *cosE = cE * (1 - 0.5 * dEA * dEA) - dEA * sE;

  } else {
    // Use Householder's third order method to guarantee performance
    // in the singular corners
    dEA = num / (denom + dEA * (0.5 * sE + one_sixth * cE * dEA));
    E = MAsign * (EA + dEA);
    *sinE = MAsign * (sE * (1 - 0.5 * dEA * dEA) + dEA * cE * (1 - dEA * dEA * one_sixth));
    *cosE = cE * (1 - 0.5 * dEA * dEA) - dEA * sE * (1 - dEA * dEA * one_sixth);
  }

  return E;
}

template <typename Scalar>
INLINE_OR_DEVICE void get_bounds(Scalar bounds[], Scalar EA_tab[], Scalar ecc) {
  const Scalar pi = 3.14159265358979323846264338327950288;
  const Scalar pi_d_12 = 3.14159265358979323846264338327950288 / 12;
  const Scalar pi_d_6 = 3.14159265358979323846264338327950288 / 6;
  const Scalar pi_d_4 = 3.14159265358979323846264338327950288 / 4;
  const Scalar pi_d_3 = 3.14159265358979323846264338327950288 / 3;
  const Scalar fivepi_d_12 = 3.14159265358979323846264338327950288 * 5. / 12;
  const Scalar pi_d_2 = 3.14159265358979323846264338327950288 / 2;
  const Scalar sevenpi_d_12 = 3.14159265358979323846264338327950288 * 7. / 12;
  const Scalar twopi_d_3 = 3.14159265358979323846264338327950288 * 2. / 3;
  const Scalar threepi_d_4 = 3.14159265358979323846264338327950288 * 3. / 4;
  const Scalar fivepi_d_6 = 3.14159265358979323846264338327950288 * 5. / 6;
  const Scalar elevenpi_d_12 = 3.14159265358979323846264338327950288 * 11. / 12;

  Scalar g2s_e = 0.2588190451025207623489 * ecc;
  Scalar g3s_e = 0.5 * ecc;
  Scalar g4s_e = 0.7071067811865475244008 * ecc;
  Scalar g5s_e = 0.8660254037844386467637 * ecc;
  Scalar g6s_e = 0.9659258262890682867497 * ecc;
  Scalar g2c_e = g6s_e;
  Scalar g3c_e = g5s_e;
  Scalar g4c_e = g4s_e;
  Scalar g5c_e = g3s_e;
  Scalar g6c_e = g2s_e;

  bounds[0] = 0;
  bounds[1] = pi_d_12 - g2s_e;
  bounds[2] = pi_d_6 - g3s_e;
  bounds[3] = pi_d_4 - g4s_e;
  bounds[4] = pi_d_3 - g5s_e;
  bounds[5] = fivepi_d_12 - g6s_e;
  bounds[6] = pi_d_2 - ecc;
  bounds[7] = sevenpi_d_12 - g6s_e;
  bounds[8] = twopi_d_3 - g5s_e;
  bounds[9] = threepi_d_4 - g4s_e;
  bounds[10] = fivepi_d_6 - g3s_e;
  bounds[11] = elevenpi_d_12 - g2s_e;
  bounds[12] = pi;

  Scalar x;

  EA_tab[1] = 1 / (1. - ecc);
  EA_tab[2] = 0;

  x = 1. / (1 - g2c_e);
  EA_tab[7] = x;
  EA_tab[8] = -0.5 * g2s_e * x * x * x;
  x = 1. / (1 - g3c_e);
  EA_tab[13] = x;
  EA_tab[14] = -0.5 * g3s_e * x * x * x;
  x = 1. / (1 - g4c_e);
  EA_tab[19] = x;
  EA_tab[20] = -0.5 * g4s_e * x * x * x;
  x = 1. / (1 - g5c_e);
  EA_tab[25] = x;
  EA_tab[26] = -0.5 * g5s_e * x * x * x;
  x = 1. / (1 - g6c_e);
  EA_tab[31] = x;
  EA_tab[32] = -0.5 * g6s_e * x * x * x;

  EA_tab[37] = 1;
  EA_tab[38] = -0.5 * ecc;

  x = 1. / (1 + g6c_e);
  EA_tab[43] = x;
  EA_tab[44] = -0.5 * g6s_e * x * x * x;
  x = 1. / (1 + g5c_e);
  EA_tab[49] = x;
  EA_tab[50] = -0.5 * g5s_e * x * x * x;
  x = 1. / (1 + g4c_e);
  EA_tab[55] = x;
  EA_tab[56] = -0.5 * g4s_e * x * x * x;
  x = 1. / (1 + g3c_e);
  EA_tab[61] = x;
  EA_tab[62] = -0.5 * g3s_e * x * x * x;
  x = 1. / (1 + g2c_e);
  EA_tab[67] = x;
  EA_tab[68] = -0.5 * g2s_e * x * x * x;

  EA_tab[73] = 1. / (1 + ecc);
  EA_tab[74] = 0;

  Scalar B0, B1, B2, idx;
  int i, k;
  for (i = 0; i < 12; i++) {
    idx = 1. / (bounds[i + 1] - bounds[i]);
    k = 6 * i;
    EA_tab[k] = i * pi_d_12;

    B0 = idx * (-EA_tab[k + 2] - idx * (EA_tab[k + 1] - idx * pi_d_12));
    B1 = idx * (-2 * EA_tab[k + 2] - idx * (EA_tab[k + 1] - EA_tab[k + 7]));
    B2 = idx * (EA_tab[k + 8] - EA_tab[k + 2]);

    EA_tab[k + 3] = B2 - 4 * B1 + 10 * B0;
    EA_tab[k + 4] = (-2 * B2 + 7 * B1 - 15 * B0) * idx;
    EA_tab[k + 5] = (B2 - 3 * B1 + 6 * B0) * idx * idx;
  }
}

template <typename Scalar>
INLINE_OR_DEVICE Scalar calcEA_fixed_ecc(const Scalar bounds[], const Scalar EA_tab[],
                                         const Scalar &M, const Scalar &ecc, Scalar *sinE,
                                         Scalar *cosE) {
  if (ecc < 1e-12 || M < 1e-12) {
    *sinE = 0;
    *cosE = 1;
    return M;
  }

  const Scalar one_sixth = 1. / 6;
  const Scalar pi = 3.14159265358979323846264338327950288;
  const Scalar pi_d_4 = 0.25 * pi;
  const Scalar pi_d_2 = 0.5 * pi;
  const Scalar threepi_d_4 = 0.75 * pi;
  const Scalar twopi = 2 * pi;

  int j, k;
  Scalar E = 0, MA, EA, sinEA, cosEA;
  Scalar dx, num, denom, dEA, dEAsq_d6;
  int MAsign = 1;
  Scalar one_over_ecc = 1e17;
  if (ecc > 1e-17) one_over_ecc = 1. / ecc;

  MA = M;
  // MA = MAmod(M);
  if (MA > pi) {
    MAsign = -1;
    MA = twopi - MA;
  }

  if (ecc < 0.78) {
    // Use the lookup table for the initial guess.
    for (j = 11; j > 0; --j)
      if (MA > bounds[j]) break;

    k = 6 * j;
    dx = MA - bounds[j];
    EA =
        EA_tab[k] + dx * (EA_tab[k + 1] +
                          dx * (EA_tab[k + 2] +
                                dx * (EA_tab[k + 3] + dx * (EA_tab[k + 4] + dx * EA_tab[k + 5]))));

    // For sinEA, since _EA in [0,pi], sinEA should always be >=0
    // (no sign ambiguity).  sqrt is much cheaper than sin.  If
    // |cos|>|sin|, compute them in reverse order, again using sqrt
    // to avoid a trig call.  Also, use trig identities, sin with a
    // low argument and the series expansion to minimize
    // computational cost.

    if (EA <= pi_d_4) {
      sinEA = shortsin(EA);
      cosEA = sqrt(1 - sinEA * sinEA);
    } else if (EA < threepi_d_4) {
      cosEA = shortsin(pi_d_2 - EA);
      sinEA = sqrt(1 - cosEA * cosEA);
    } else {
      sinEA = shortsin(pi - EA);
      cosEA = -sqrt(1 - sinEA * sinEA);
    }

    num = (MA - EA) * one_over_ecc + sinEA;
    denom = one_over_ecc - cosEA;

    // Second order approximation.
    dEA = num * denom / (denom * denom + 0.5 * sinEA * num);

    // Apply our correction to EA, sinEA, and cosEA using
    // series.  Go to second order, since that was our level of
    // approximation above and will get us to basically machine
    // precision for eccentricities below 0.78.

    E = MAsign * (EA + dEA);
    *sinE = MAsign * (sinEA * (1 - 0.5 * dEA * dEA) + dEA * cosEA);
    *cosE = cosEA * (1 - 0.5 * dEA * dEA) - dEA * sinEA;
  } else {
    // Higher eccentricities will require a third-order correction to
    // achieve machine precision for all values of the eccentric
    // anomaly.  In the singular corner, they also use a series
    // expansion rather than the piecewise polynomial fit.

    if (2 * MA + (1 - ecc) > 0.2) {
      // Use the lookup table for the initial guess as long as we
      // are not in the singular corner.
      for (j = 11; j > 0; --j)
        if (MA > bounds[j]) break;

      k = 6 * j;
      dx = MA - bounds[j];
      EA = EA_tab[k] +
           dx * (EA_tab[k + 1] +
                 dx * (EA_tab[k + 2] +
                       dx * (EA_tab[k + 3] + dx * (EA_tab[k + 4] + dx * EA_tab[k + 5]))));
    } else {
      // Use the series expansions in the singular corner.
      EA = EAstart(MA, ecc);
    }

    if (EA <= pi_d_4) {
      sinEA = shortsin(EA);
      cosEA = sqrt(1 - sinEA * sinEA);
    } else if (EA < threepi_d_4) {
      cosEA = shortsin(pi_d_2 - EA);
      sinEA = sqrt(1 - cosEA * cosEA);
    } else {
      sinEA = shortsin(pi - EA);
      cosEA = -sqrt(1 - sinEA * sinEA);
    }

    num = (MA - EA) * one_over_ecc + sinEA;
    denom = one_over_ecc - cosEA;

    if (MA > 0.4) {
      dEA = num * denom / (denom * denom + 0.5 * sinEA * num);

    } else {
      dEA = num * (denom * denom + 0.5 * num * sinEA);
      dEA /= denom * denom * denom + num * (denom * sinEA + one_sixth * num * cosEA);
    }

    dEAsq_d6 = dEA * dEA * one_sixth;

    // Apply our correction to EA, sinEA, and cosEA using
    // series.  Go to third order, since that was our level of
    // approximation above and will get us to basically machine
    // precision at the higher eccentricities.
    E = MAsign * (EA + dEA);
    *sinE = MAsign * (sinEA * (1 - 3 * dEAsq_d6) + dEA * (1 - dEAsq_d6) * cosEA);
    *cosE = cosEA * (1 - 3 * dEAsq_d6) - dEA * (1 - dEAsq_d6) * sinEA;
  }
  return E;
}

}  // namespace brandt

template <typename T>
struct Brandt : public Solver<T> {
  INLINE_OR_DEVICE T compute_eccentric_anomaly(const T &mean_anomaly,
                                               const T &eccentricity) override {
    T sinE, cosE;
    return brandt::calcEA(mean_anomaly, eccentricity, &sinE, &cosE);
  }

  INLINE_OR_DEVICE void evaluate(const T &mean_anomaly, const T &eccentricity, T *sin_ecc_anomaly,
                                 T *cos_ecc_anomaly) override {
    brandt::calcEA(mean_anomaly, eccentricity, sin_ecc_anomaly, cos_ecc_anomaly);
  }
};

template <typename T>
struct BrandtFixed : public Solver<T> {
  T EA_tab[6 * 13], bounds[13];

  INLINE_OR_DEVICE void precompute_for_eccentricity(const T &eccentricity) override {
    this->eccentricity_ = eccentricity;
    brandt::get_bounds(bounds, EA_tab, eccentricity);
  }

  INLINE_OR_DEVICE T compute_eccentric_anomaly(const T &mean_anomaly) override {
    T sinE, cosE;
    return brandt::calcEA_fixed_ecc(bounds, EA_tab, mean_anomaly, this->eccentricity_, &sinE,
                                    &cosE);
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
    brandt::calcEA_fixed_ecc(bounds, EA_tab, mean_anomaly, this->eccentricity_, sin_ecc_anomaly,
                             cos_ecc_anomaly);
  }
};

}  // namespace kepler_benchmarks

#endif