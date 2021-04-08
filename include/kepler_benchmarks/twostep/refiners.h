#ifndef _KEPLER_BENCHMARKS_TWOSTEP_REFINERS_H_
#define _KEPLER_BENCHMARKS_TWOSTEP_REFINERS_H_

#include "kepler_benchmarks/cuda_support.h"
#include "kepler_benchmarks/math_utils.h"

namespace kepler_benchmarks {
namespace twostep {

template <typename T>
struct Refiner {
  virtual INLINE_OR_DEVICE T refine_estimate(const T& reduced_mean_anomaly, const T& eccentricity,
                                             const T& eccentric_anomaly, T* sin_eccentric_anomaly,
                                             T* cos_eccentric_anomaly) = 0;
};

// Nijenhuis (1991)
// http://adsabs.harvard.edu/abs/1991CeMDA..51..319N
namespace nijenhuis {

// Calculates x - sin(x) and 1 - cos(x) to 20 significant digits for x in [0, pi)
template <typename T>
inline void sin_cos_reduc(T x, T* SnReduc, T* CsReduc) {
  const T s[] = {1.0 / 6,   1.0 / 20,  1.0 / 42,  1.0 / 72,  1.0 / 110,
                 1.0 / 156, 1.0 / 210, 1.0 / 272, 1.0 / 342, 1.0 / 420};
  const T c[] = {0.5,       1.0 / 12,  1.0 / 30,  1.0 / 56,  1.0 / 90,
                 1.0 / 132, 1.0 / 182, 1.0 / 240, 1.0 / 306, 1.0 / 380};

  bool bigg = x > M_PI_2;
  T u = (bigg) ? M_PI - x : x;
  bool big = u > M_PI_2;
  T v = (big) ? M_PI_2 - u : u;
  T w = v * v;

  T ss = T(1);
  T cc = T(1);
  for (int i = 9; i >= 1; --i) {
    ss = 1 - w * s[i] * ss;
    cc = 1 - w * c[i] * cc;
  }
  ss *= v * w * s[0];
  cc *= w * c[0];

  if (big) {
    *SnReduc = u - 1 + cc;
    *CsReduc = 1 - M_PI_2 + u + ss;
  } else {
    *SnReduc = ss;
    *CsReduc = cc;
  }
  if (bigg) {
    *SnReduc = 2 * x - M_PI + *SnReduc;
    *CsReduc = 2 - *CsReduc;
  }
}

}  // namespace nijenhuis

template <typename T>
struct NijenhuisRefiner : public Starter<T> {
  INLINE_OR_DEVICE T refine_estimate(const T& reduced_mean_anomaly, const T& eccentricity,
                                     const T& eccentric_anomaly, T* sin_eccentric_anomaly,
                                     T* cos_eccentric_anomaly) override {
    T sE, cE, E = eccentric_anomaly;
    nijenhuis::sin_cos_reduc(E, &sE, &cE);
    // T sE = E - sin(E);
    // T cE = 1 - cos(E);

    const T ome = 1 - eccentricity;
    const T f_0 = eccentricity * sE + E * ome - M;
    const T f_1 = eccentricity * cE + ome;
    const T f_2 = eccentricity * (E - sE);
    const T f_3 = 1 - f_1;
    const T d_3 = -f_0 / (f_1 - 0.5 * f_0 * f_2 / f_1);
    const T d_4 = -f_0 / (f_1 + 0.5 * d_3 * f_2 + (d_3 * d_3) * f_3 / 6);
    const T d_42 = d_4 * d_4;
    const T dE = -f_0 / (f_1 + 0.5 * d_4 * f_2 + d_4 * d_4 * f_3 / 6 - d_42 * d_4 * f_2 / 24);
    E += dE;

    *sin_eccentric_anomaly = sin(E);
    *cos_eccentric_anomaly = cos(E);
    return E;
  }
};

}  // namespace twostep
}  // namespace kepler_benchmarks

#endif