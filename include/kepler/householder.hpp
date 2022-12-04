#ifndef KEPLER_HOUSEHOLDER_HPP
#define KEPLER_HOUSEHOLDER_HPP

#include "./constants.hpp"
#include "./math.hpp"
#include "./simd.hpp"

namespace kepler {
namespace householder {
namespace detail {

// https://aherrmann.github.io/programming/2016/02/28/unpacking-tuples-in-cpp14/
template <typename Tuple, size_t... Is>
constexpr auto tuple_init_impl(Tuple t, std::index_sequence<Is...>) {
  return std::make_tuple(std::get<Is>(t)...);
}

template <typename Tuple>
constexpr auto tuple_init(Tuple t) {
  return tuple_init_impl(t, std::make_index_sequence<std::tuple_size<Tuple>{} - 1>{});
}

template <typename T>
inline T householder(const T& f0, const T& f1) {
  return -f0 / f1;
}

template <typename T>
inline T householder_impl(const T& f0, const T& f1, const std::tuple<T>&) {
  auto d = householder(f0, f1);
  return -f0 / math::horner(d, f1);
}

template <typename T, typename Args>
inline T householder_impl(const T& f0, const T& f1, const Args& args) {
  auto d = householder_impl(f0, f1, tuple_init(args));
  return -f0 / math::horner_impl(d, f1, args);
}

template <typename T, typename... Args>
inline T householder(const T& f0, const T& f1, Args... args) {
  return householder_impl(f0, f1, std::make_tuple(args...));
}

template <typename T>
struct householder_state {
  T f0;
  T sin;
  T cos;
};

}  // namespace detail

template <int order>
struct householder {
  static_assert(order > 0, "order must be positive");
  static_assert(order <= 7, "order cannot be larger than 7");

  template <typename T>
  static inline detail::householder_state<T> init(const T& eccentricity, const T& mean_anomaly,
                                                  T& eccentric_anomaly) {
    auto sincos = math::sincos(eccentric_anomaly);
    auto f0 = eccentric_anomaly - math::fma(eccentricity, sincos.first, mean_anomaly);
    return detail::householder_state<T>{f0, sincos.first, sincos.second};
  }

  template <typename A, typename T>
  static inline detail::householder_state<xs::batch<T, A>> init(
      const T& eccentricity, const xs::batch<T, A>& mean_anomaly,
      xs::batch<T, A>& eccentric_anomaly) {
    using B = xs::batch<T, A>;
    auto sincos = math::sincos(eccentric_anomaly);
    auto f0 = eccentric_anomaly - xs::fma(B(eccentricity), sincos.first, mean_anomaly);
    return detail::householder_state<B>{f0, sincos.first, sincos.second};
  }

  template <typename T, typename B>
  static inline B step(const detail::householder_state<B>& state, const T& eccentricity) {
    auto f0 = state.f0;
    auto s = state.sin;
    auto c = state.cos;

    auto f1 = math::fnma<B>(B(eccentricity), c, B(T(1.)));
    KEPLER_IF_CONSTEXPR(order == 1) { return detail::householder(f0, f1); }

    auto f2 = eccentricity * s;
    KEPLER_IF_CONSTEXPR(order == 2) {
      return detail::householder(f0, f1, constants::hh2<T>() * f2);
    }

    auto f3 = T(1.) - f1;
    KEPLER_IF_CONSTEXPR(order == 3) {
      return detail::householder(f0, f1, constants::hh2<T>() * f2, constants::hh3<T>() * f3);
    }

    auto f4 = -f2;
    KEPLER_IF_CONSTEXPR(order == 4) {
      return detail::householder(f0, f1, constants::hh2<T>() * f2, constants::hh3<T>() * f3,
                                 constants::hh4<T>() * f4);
    }

    auto f5 = -f3;
    KEPLER_IF_CONSTEXPR(order == 5) {
      return detail::householder(f0, f1, constants::hh2<T>() * f2, constants::hh3<T>() * f3,
                                 constants::hh4<T>() * f4, constants::hh5<T>() * f5);
    }

    auto f6 = -f4;
    KEPLER_IF_CONSTEXPR(order == 6) {
      return detail::householder(f0, f1, constants::hh2<T>() * f2, constants::hh3<T>() * f3,
                                 constants::hh4<T>() * f4, constants::hh5<T>() * f5,
                                 constants::hh6<T>() * f6);
    }

    auto f7 = -f5;
    return detail::householder(f0, f1, constants::hh2<T>() * f2, constants::hh3<T>() * f3,
                               constants::hh4<T>() * f4, constants::hh5<T>() * f5,
                               constants::hh6<T>() * f6, constants::hh7<T>() * f7);
  }
};

}  // namespace householder
}  // namespace kepler

#endif