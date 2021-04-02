#ifndef _KEPLER_BENCHMARKS_SOLVER_H_
#define _KEPLER_BENCHMARKS_SOLVER_H_

#include "cuda_support.h"

namespace kepler_benchmarks {

template <typename T>
struct Solver {
  typedef T Scalar;

  T eccentricity_;

  Solver() : eccentricity_(T(0)){};

  virtual INLINE_OR_DEVICE void precompute_for_eccentricity(const T& eccentricity) {
    eccentricity_ = eccentricity;
  }

  virtual INLINE_OR_DEVICE T compute_eccentric_anomaly(const T& mean_anomaly,
                                                       const T& eccentricity) = 0;

  virtual INLINE_OR_DEVICE T compute_eccentric_anomaly(const T& mean_anomaly) {
    return this->compute_eccentric_anomaly(mean_anomaly, eccentricity_);
  }

  virtual INLINE_OR_DEVICE void evaluate(const T& mean_anomaly, const T& eccentricity,
                                         T* sin_ecc_anomaly, T* cos_ecc_anomaly) {
    T E = this->compute_eccentric_anomaly(mean_anomaly, eccentricity);
    *sin_ecc_anomaly = sin(E);
    *cos_ecc_anomaly = cos(E);
  };

  virtual INLINE_OR_DEVICE void evaluate(const T& mean_anomaly, T* sin_ecc_anomaly,
                                         T* cos_ecc_anomaly) {
    return this->evaluate(mean_anomaly, eccentricity_, sin_ecc_anomaly, cos_ecc_anomaly);
  }
};

template <typename T>
struct DummyBaselineSolver : public Solver<T> {
  virtual INLINE_OR_DEVICE T compute_eccentric_anomaly(const T& mean_anomaly,
                                                       const T& eccentricity) override {
    (void)(eccentricity);
    return mean_anomaly;
  };
};

}  // namespace kepler_benchmarks

#endif