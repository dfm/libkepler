#ifndef KEPLER_HOUSEHOLDER_HPP
#define KEPLER_HOUSEHOLDER_HPP

#include <cstdint>

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

template <typename Tuple, typename Is = std::make_index_sequence<std::tuple_size<Tuple>{} - 1>>
constexpr auto tuple_init(Tuple t) {
  return tuple_init_impl(t, Is{});
}

template <std::uint32_t n>
struct factorial {
  constexpr static std::uint32_t value = n * factorial<n - 1>::value;
};

template <>
struct factorial<0> {
  constexpr static std::uint32_t value = 1;
};

template <typename T, typename Tuple, size_t... Is>
inline T horner_packed_impl(const T& x, Tuple t, std::index_sequence<Is...>) {
  return math::horner(x, std::get<Is>(t)...);
}

template <typename T>
inline T horner_packed(const T&, const std::tuple<T>& t) {
  return std::get<0>(t);
}

template <typename T, typename Tuple>
inline T horner_packed(const T& x, Tuple t) {
  return horner_packed_impl(x, t, std::make_index_sequence<std::tuple_size<Tuple>{}>{});
}

template <typename T>
inline T householder_packed(const T& f0, const std::tuple<T>& args) {
  return -f0 / std::get<0>(args);
}

template <typename T, typename Args>
inline T householder_packed(const T& f0, const Args& args) {
  auto d = householder_packed(f0, tuple_init(args));
  return -f0 / horner_packed(d, args);
}

template <typename T, typename Tuple, size_t... Is>
inline T householder_impl(const T& f0, Tuple t, std::index_sequence<Is...>) {
  return householder_packed(f0, std::make_tuple(std::get<Is>(t) / factorial<Is + 1>::value...));
}

template <typename T, typename... Args>
inline T householder(const T& f0, Args... args) {
  return householder_impl(f0, std::make_tuple(args...), std::index_sequence_for<Args...>{});
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
      return detail::householder(f0, f1, f2);
      // return detail::householder(f0, f1, constants::hh2<T>() * f2);
    }

    auto f3 = T(1.) - f1;
    KEPLER_IF_CONSTEXPR(order == 3) {
      return detail::householder(f0, f1, f2, f3);
      // return detail::householder(f0, f1, constants::hh2<T>() * f2, constants::hh3<T>() * f3);
    }

    auto f4 = -f2;
    KEPLER_IF_CONSTEXPR(order == 4) {
      return detail::householder(f0, f1, f2, f3, f4);
      // return detail::householder(f0, f1, constants::hh2<T>() * f2, constants::hh3<T>() * f3,
      //                            constants::hh4<T>() * f4);
    }

    auto f5 = -f3;
    KEPLER_IF_CONSTEXPR(order == 5) {
      return detail::householder(f0, f1, f2, f3, f4, f5);
      // return detail::householder(f0, f1, constants::hh2<T>() * f2, constants::hh3<T>() * f3,
      //                            constants::hh4<T>() * f4, constants::hh5<T>() * f5);
    }

    auto f6 = -f4;
    KEPLER_IF_CONSTEXPR(order == 6) {
      return detail::householder(f0, f1, f2, f3, f4, f5, f6);
      // return detail::householder(f0, f1, constants::hh2<T>() * f2, constants::hh3<T>() * f3,
      //                            constants::hh4<T>() * f4, constants::hh5<T>() * f5,
      //                            constants::hh6<T>() * f6);
    }

    auto f7 = -f5;
    return detail::householder(f0, f1, f2, f3, f4, f5, f6, f7);
    // return detail::householder(f0, f1, constants::hh2<T>() * f2, constants::hh3<T>() * f3,
    //                            constants::hh4<T>() * f4, constants::hh5<T>() * f5,
    //                            constants::hh6<T>() * f6, constants::hh7<T>() * f7);
  }
};

}  // namespace householder
}  // namespace kepler

#endif