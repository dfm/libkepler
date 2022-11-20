#ifndef KEPLER_REFINERS_HPP
#define KEPLER_REFINERS_HPP

#include <cmath>

#include "./constants.hpp"
#include "./householder.hpp"
#include "./starters.hpp"

namespace kepler {
namespace refiners {

namespace detail {
template <typename T>
inline T default_tolerance() {
  return T(1e-12);
}

template <>
inline float default_tolerance() {
  return 1e-6f;
}
}  // namespace detail

template <int order, typename T>
struct iterative {
  typedef T value_type;
  typedef starters::basic<T> default_starter;

  int max_iterations;
  T tolerance;
  iterative() : max_iterations(30), tolerance(detail::default_tolerance<T>()) {}
  iterative(T tolerance) : max_iterations(30), tolerance(tolerance) {}
  iterative(int max_iterations, T tolerance)
      : max_iterations(max_iterations), tolerance(tolerance) {}

  inline T refine(const T& eccentricity, const T& mean_anomaly,
                  const T& initial_eccentric_anomaly) const {
    T eccentric_anomaly = initial_eccentric_anomaly;
    for (int i = 0; i < max_iterations; ++i) {
      auto delta = ::kepler::detail::householder<order>::step(eccentricity, mean_anomaly,
                                                              eccentric_anomaly);
      if (std::abs(delta) < tolerance) break;
    }
    return eccentric_anomaly;
  }
};

template <int order, typename T>
struct non_iterative {
  typedef T value_type;
  typedef starters::basic<T> default_starter;
  inline T refine(const T& eccentricity, const T& mean_anomaly,
                  const T& initial_eccentric_anomaly) const {
    T eccentric_anomaly = initial_eccentric_anomaly;
    ::kepler::detail::householder<order>::step(eccentricity, mean_anomaly, eccentric_anomaly);
    return eccentric_anomaly;
  }
};

}  // namespace refiners
}  // namespace kepler

#endif
