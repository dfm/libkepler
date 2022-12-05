#ifndef KEPLER_REFINERS_HPP
#define KEPLER_REFINERS_HPP

#include <cmath>

#include "./constants.hpp"
#include "./householder.hpp"
#include "./simd.hpp"
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

template <typename T>
struct noop {
  typedef T value_type;

  inline T refine(const T&, const T&, const T& initial_eccentric_anomaly) const {
    return initial_eccentric_anomaly;
  }

  template <typename A>
  inline xs::batch<T, A> refine(const T&, const xs::batch<T, A>&,
                                const xs::batch<T, A>& initial_eccentric_anomaly) const {
    return initial_eccentric_anomaly;
  }
};

template <int order, typename T>
struct iterative {
  typedef T value_type;
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
      auto state = householder::init(eccentricity, mean_anomaly, eccentric_anomaly);
      if (std::abs(state.f0) < tolerance) break;
      eccentric_anomaly += householder::step<order>(state);
    }
    return eccentric_anomaly;
  }

  template <typename A>
  inline xs::batch<T, A> refine(const T& eccentricity, const xs::batch<T, A>& mean_anomaly,
                                const xs::batch<T, A>& initial_eccentric_anomaly) const {
    using B = xs::batch<T, A>;
    B eccentric_anomaly = initial_eccentric_anomaly;
    typename B::batch_bool_type converged(false);
    for (int i = 0; i < max_iterations; ++i) {
      auto state = householder::init(eccentricity, mean_anomaly, eccentric_anomaly);
      converged = converged | (xs::abs(state.f0) < B(tolerance));
      if (xs::all(converged)) break;
      auto delta = householder::step<order>(state);
      eccentric_anomaly = xs::select(converged, eccentric_anomaly, eccentric_anomaly + delta);
    }
    return eccentric_anomaly;
  }
};

template <int order, typename T>
struct non_iterative {
  typedef T value_type;

  template <typename B>
  inline B refine(const T& eccentricity, const B& mean_anomaly,
                  const B& initial_eccentric_anomaly) const {
    B eccentric_anomaly = initial_eccentric_anomaly;
    auto state = householder::init(eccentricity, mean_anomaly, eccentric_anomaly);
    eccentric_anomaly += householder::step<order>(state);
    return eccentric_anomaly;
  }
};

}  // namespace refiners
}  // namespace kepler

#endif
