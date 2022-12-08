#ifndef KEPLER_REFINERS_HPP
#define KEPLER_REFINERS_HPP

#include <cmath>

#include "constants.hpp"
#include "householder.hpp"
#include "simd.hpp"
#include "starters.hpp"

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

template <typename T>
struct _refiner {
  typedef T value_type;
};
}  // namespace detail

template <typename T>
struct noop : detail::_refiner<T> {
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
struct iterative : detail::_refiner<T> {
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

namespace detail {

template <int order, typename T, typename B>
static inline B _non_iterative_step(const T& eccentricity, const B& mean_anomaly,
                                    const B& initial_eccentric_anomaly) {
  auto state = householder::init(eccentricity, mean_anomaly, initial_eccentric_anomaly);
  return initial_eccentric_anomaly + householder::step<order>(state);
}

template <size_t level>
struct _non_iterative {
  template <int order, typename T, typename B>
  static inline B refine(const T& eccentricity, const B& mean_anomaly,
                         const B& initial_eccentric_anomaly) {
    B eccentric_anomaly =
        _non_iterative_step<order, T, B>(eccentricity, mean_anomaly, initial_eccentric_anomaly);
    return _non_iterative<level - 1>::template refine<order, T, B>(eccentricity, mean_anomaly,
                                                                   eccentric_anomaly);
  }
};

template <>
struct _non_iterative<1> {
  template <int order, typename T, typename B>
  static inline B refine(const T& eccentricity, const B& mean_anomaly,
                         const B& initial_eccentric_anomaly) {
    return _non_iterative_step<order, T, B>(eccentricity, mean_anomaly, initial_eccentric_anomaly);
  }
};

}  // namespace detail

template <int order, typename T, size_t num = 1>
struct non_iterative : detail::_refiner<T> {
  template <typename B>
  inline B refine(const T& eccentricity, const B& mean_anomaly,
                  const B& initial_eccentric_anomaly) const {
    return detail::_non_iterative<num>::template refine<order, T, B>(eccentricity, mean_anomaly,
                                                                     initial_eccentric_anomaly);
  }
};

template <typename T>
struct brandt : detail::_refiner<T> {
  inline T refine(const T& eccentricity, const T& mean_anomaly,
                  const T& initial_eccentric_anomaly) const {
    if (eccentricity < T(0.78) || mean_anomaly > T(0.4)) {
      return detail::_non_iterative_step<2>(eccentricity, mean_anomaly, initial_eccentric_anomaly);
    } else {
      return detail::_non_iterative_step<3>(eccentricity, mean_anomaly, initial_eccentric_anomaly);
    }
  }

  template <typename A>
  inline xs::batch<T, A> refine(const T& eccentricity, const xs::batch<T, A>& mean_anomaly,
                                const xs::batch<T, A>& initial_eccentric_anomaly) const {
    using B = xs::batch<T, A>;
    auto flag = typename B::batch_bool_type(eccentricity < T(0.78)) | (mean_anomaly > B(T(0.4)));
    if (xs::all(flag)) {
      return detail::_non_iterative_step<2>(eccentricity, mean_anomaly, initial_eccentric_anomaly);
    } else if (xs::none(flag)) {
      return detail::_non_iterative_step<3>(eccentricity, mean_anomaly, initial_eccentric_anomaly);
    } else {
      auto fast =
          detail::_non_iterative_step<2>(eccentricity, mean_anomaly, initial_eccentric_anomaly);
      auto slow =
          detail::_non_iterative_step<3>(eccentricity, mean_anomaly, initial_eccentric_anomaly);
      return xs::select(flag, fast, slow);
    }
  }
};

template <typename R>
struct refine_with_eccentricity : detail::_refiner<typename R::value_type> {
  using T = typename R::value_type;

  template <typename B>
  static inline B refine(const R& refiner, const T& eccentricity, const B& mean_anomaly,
                         const B& initial_eccentric_anomaly, B* sin_eccentric_anomaly,
                         B* cos_eccentric_anomaly) {
    auto ecc_anom = refiner.refine(eccentricity, mean_anomaly, initial_eccentric_anomaly);
    auto sincos = math::sincos(ecc_anom);
    *sin_eccentric_anomaly = sincos.first;
    *cos_eccentric_anomaly = sincos.second;
    return ecc_anom;
  };
};

template <typename T>
struct refine_with_eccentricity<brandt<T>> : detail::_refiner<T> {
  static inline T refine(const brandt<T>&, const T& eccentricity, const T& mean_anomaly,
                         const T& initial_eccentric_anomaly, T* sin_eccentric_anomaly,
                         T* cos_eccentric_anomaly) {
    if (eccentricity < detail::default_tolerance<T>()) {
      auto sincos = math::sincos(initial_eccentric_anomaly);
      *sin_eccentric_anomaly = sincos.first;
      *cos_eccentric_anomaly = sincos.second;
      return initial_eccentric_anomaly;
    } else if (eccentricity < T(0.78) || mean_anomaly > T(0.4)) {
      auto state = householder::init(eccentricity, mean_anomaly, initial_eccentric_anomaly);
      auto delta = householder::step<2>(state);
      auto factor = math::fma(T(-0.5) * delta, delta, T(1.));
      *sin_eccentric_anomaly =
          math::fma(factor, state.ecc_sin, delta * state.ecc_cos) / eccentricity;
      *cos_eccentric_anomaly =
          math::fnma(delta, state.ecc_sin, factor * state.ecc_cos) / eccentricity;
      return initial_eccentric_anomaly + delta;
    } else {
      auto state = householder::init(eccentricity, mean_anomaly, initial_eccentric_anomaly);
      auto delta = householder::step<3>(state);
      auto factor = constants::sixth<T>() * delta * delta;
      auto factor1 = math::fnma(T(3.), factor, T(1.));
      auto factor2 = math::fnma(delta, factor, delta);
      *sin_eccentric_anomaly =
          math::fma(factor1, state.ecc_sin, factor2 * state.ecc_cos) / eccentricity;
      *cos_eccentric_anomaly =
          math::fnma(factor2, state.ecc_sin, factor1 * state.ecc_cos) / eccentricity;
      return initial_eccentric_anomaly + delta;
    }
  }

  // TODO(dfm): Update the batch version to update sin and cos properly too!
  template <typename A>
  static inline xs::batch<T, A> refine(const brandt<T>& refiner, const T& eccentricity,
                                       const xs::batch<T, A>& mean_anomaly,
                                       const xs::batch<T, A>& initial_eccentric_anomaly,
                                       xs::batch<T, A>* sin_eccentric_anomaly,
                                       xs::batch<T, A>* cos_eccentric_anomaly) {
    auto ecc_anom = refiner.refine(eccentricity, mean_anomaly, initial_eccentric_anomaly);
    auto sincos = math::sincos(ecc_anom);
    *sin_eccentric_anomaly = sincos.first;
    *cos_eccentric_anomaly = sincos.second;
    return ecc_anom;
  }
};

}  // namespace refiners
}  // namespace kepler

#endif
