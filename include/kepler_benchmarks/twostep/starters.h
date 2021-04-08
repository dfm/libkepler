#ifndef _KEPLER_BENCHMARKS_TWOSTEP_STARTERS_H_
#define _KEPLER_BENCHMARKS_TWOSTEP_STARTERS_H_

#include "kepler_benchmarks/cuda_support.h"
#include "kepler_benchmarks/math_utils.h"

namespace kepler_benchmarks {
namespace twostep {

template <typename T>
struct Starter {
  T eccentricity_;

  Starter() : eccentricity_(T(0)){};

  virtual INLINE_OR_DEVICE void precompute_for_eccentricity(const T& eccentricity) {
    eccentricity_ = eccentricity;
  }

  virtual INLINE_OR_DEVICE T get_starter_reduced(const T& mean_anomaly) = 0;
};

// Markley (1995)
// http://adsabs.harvard.edu/abs/1995CeMDA..63..101M
template <typename T>
struct MarkleyStarter : public Starter<T> {
  INLINE_OR_DEVICE T get_starter_reduced(const T& reduced_mean_anomaly) override {
    // M must be in the range [0, pi)
    const T FACTOR1 = 3 * M_PI / (M_PI - 6 / M_PI);
    const T FACTOR2 = 1.6 / (M_PI - 6 / M_PI);

    const T M = reduced_mean_anomaly;
    const T ecc = eccentricity_;
    const T ome = 1 - ecc;

    const T M2 = M * M;
    const T alpha = FACTOR1 + FACTOR2 * (M_PI - M) / (1 + ecc);
    const T d = 3 * ome + alpha * ecc;
    const T alphad = alpha * d;
    const T r = (3 * alphad * (d - ome) + M2) * M;
    const T q = 2 * alphad * ome - M2;
    const T q2 = q * q;

    T w = cbrt(std::abs(r) + sqrt(q2 * q + r * r));
    w *= w;

    return (2 * r * w / (w * w + w * q + q2) + M) / d;
  }
};

}  // namespace twostep
}  // namespace kepler_benchmarks

#endif