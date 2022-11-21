#ifndef KEPLER_HOUSEHOLDER_HPP
#define KEPLER_HOUSEHOLDER_HPP

#include "./constants.hpp"
#include "./simd.hpp"

namespace kepler {
namespace detail {

template <typename T>
inline T householder_eval(const T& f0, const T& f1) {
  return -f0 / f1;
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2) {
  auto d = householder_eval(f0, f1);
  return -f0 / (f1 + d * constants::hh2<T>() * f2);
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2, const T& f3) {
  auto d = householder_eval(f0, f1, f2);
  return -f0 / (f1 + d * (constants::hh2<T>() * f2 + d * constants::hh3<T>() * f3));
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2, const T& f3, const T& f4) {
  auto d = householder_eval(f0, f1, f2, f3);
  return -f0 / (f1 + d * (constants::hh2<T>() * f2 +
                          d * (constants::hh3<T>() * f3 + d * constants::hh4<T>() * f4)));
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2, const T& f3, const T& f4,
                          const T& f5) {
  auto d = householder_eval(f0, f1, f2, f3);
  return -f0 / (f1 + d * (constants::hh2<T>() * f2 +
                          d * (constants::hh3<T>() * f3 +
                               d * (constants::hh4<T>() * f4 + d * constants::hh5<T>() * f5))));
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2, const T& f3, const T& f4,
                          const T& f5, const T& f6) {
  auto d = householder_eval(f0, f1, f2, f3);
  return -f0 /
         (f1 + d * (constants::hh2<T>() * f2 +
                    d * (constants::hh3<T>() * f3 +
                         d * (constants::hh4<T>() * f4 +
                              d * (constants::hh5<T>() * f5 + d * constants::hh6<T>() * f6)))));
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2, const T& f3, const T& f4,
                          const T& f5, const T& f6, const T& f7) {
  auto d = householder_eval(f0, f1, f2, f3);
  return -f0 /
         (f1 +
          d * (constants::hh2<T>() * f2 +
               d * (constants::hh3<T>() * f3 +
                    d * (constants::hh4<T>() * f4 +
                         d * (constants::hh5<T>() * f5 +
                              d * (constants::hh6<T>() * f6 + d * constants::hh7<T>() * f7))))));
}

template <typename T>
struct householder_state {
  T f0;
  T sin;
  T cos;
};

template <int order>
struct householder {
  static_assert(order > 0, "order must be positive");
  static_assert(order <= 7, "order cannot be larger than 7");

  template <typename T>
  static inline householder_state<T> init(const T& eccentricity, const T& mean_anomaly,
                                          T& eccentric_anomaly) {
    auto s = std::sin(eccentric_anomaly);
    auto c = std::cos(eccentric_anomaly);
    auto f0 = eccentric_anomaly - eccentricity * s - mean_anomaly;
    return householder_state<T>{f0, s, c};
  }

  template <typename A, typename T>
  static inline householder_state<xs::batch<T, A>> init(const T& eccentricity,
                                                        const xs::batch<T, A>& mean_anomaly,
                                                        xs::batch<T, A>& eccentric_anomaly) {
    using B = xs::batch<T, A>;
    auto sincos = xs::sincos(eccentric_anomaly);
    auto f0 = eccentric_anomaly - xs::fma(B(eccentricity), sincos.first, mean_anomaly);
    return householder_state<B>{f0, sincos.first, sincos.second};
  }

  template <typename T>
  static inline T step(const householder_state<T>& state, const T& eccentricity) {
    auto f0 = state.f0;
    auto s = state.sin;
    auto c = state.cos;

    auto f1 = T(1.) - eccentricity * c;
    KEPLER_IF_CONSTEXPR(order == 1) { return householder_eval(f0, f1); }

    auto f2 = eccentricity * s;
    KEPLER_IF_CONSTEXPR(order == 2) { return householder_eval(f0, f1, f2); }

    auto f3 = T(1.) - f1;
    KEPLER_IF_CONSTEXPR(order == 3) { return householder_eval(f0, f1, f2, f3); }

    auto f4 = -f2;
    KEPLER_IF_CONSTEXPR(order == 4) { return householder_eval(f0, f1, f2, f3, f4); }

    auto f5 = -f3;
    KEPLER_IF_CONSTEXPR(order == 5) { return householder_eval(f0, f1, f2, f3, f4, f5); }

    auto f6 = -f4;
    KEPLER_IF_CONSTEXPR(order == 5) { return householder_eval(f0, f1, f2, f3, f4, f5, f6); }

    auto f7 = -f5;
    return householder_eval(f0, f1, f2, f3, f4, f5, f6, f7);
  }

  template <typename A, typename T>
  static inline xs::batch<T, A> step(const householder_state<xs::batch<T, A>>& state,
                                     const T& eccentricity) {
    using B = xs::batch<T, A>;
    auto f0 = state.f0;
    auto s = state.sin;
    auto c = state.cos;

    auto f1 = xs::fnma(B(eccentricity), c, B(T(1.)));
    KEPLER_IF_CONSTEXPR(order == 1) { return householder_eval(f0, f1); }

    auto f2 = eccentricity * s;
    KEPLER_IF_CONSTEXPR(order == 2) { return householder_eval(f0, f1, f2); }

    auto f3 = B(T(1.)) - f1;
    KEPLER_IF_CONSTEXPR(order == 3) { return householder_eval(f0, f1, f2, f3); }

    auto f4 = -f2;
    KEPLER_IF_CONSTEXPR(order == 4) { return householder_eval(f0, f1, f2, f3, f4); }

    auto f5 = -f3;
    KEPLER_IF_CONSTEXPR(order == 5) { return householder_eval(f0, f1, f2, f3, f4, f5); }

    auto f6 = -f4;
    KEPLER_IF_CONSTEXPR(order == 5) { return householder_eval(f0, f1, f2, f3, f4, f5, f6); }

    auto f7 = -f5;
    return householder_eval(f0, f1, f2, f3, f4, f5, f6, f7);
  }
};

}  // namespace detail
}  // namespace kepler

#endif