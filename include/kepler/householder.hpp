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

template <std::uint64_t n>
struct factorial {
  constexpr static std::uint64_t value = n * factorial<n - 1>::value;
};

template <>
struct factorial<0> {
  constexpr static std::uint64_t value = 1;
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
struct state {
  T f0;
  T ecc_sin;
  T ecc_cos;
};
template <bool is_even>
struct evaluate_impl {
  template <typename T>
  static inline T value(const state<T>& state) {
    return state.ecc_sin;
  }
};

template <>
struct evaluate_impl<false> {
  template <typename T>
  static inline T value(const state<T>& state) {
    return state.ecc_cos;
  }
};

template <size_t order>
struct evaluate {
  static_assert(order > 0, "order must be larger than 0");
  template <typename T>
  static inline T get(const state<T>& s) {
    constexpr int sign = order % 4 < 2 ? -1 : 1;
    return T(sign) * evaluate_impl<order % 2 == 0>::value(s);
  }
};

template <>
struct evaluate<1> {
  template <typename T>
  static inline T get(const state<T>& s) {
    return T(1.) - s.ecc_cos;
  }
};

template <>
struct evaluate<2> {
  template <typename T>
  static inline T get(const state<T>& s) {
    return s.ecc_sin;
  }
};

template <>
struct evaluate<3> {
  template <typename T>
  static inline T get(const state<T>& s) {
    return s.ecc_cos;
  }
};

template <typename T, typename Seq>
struct expander;

template <typename T, std::size_t... Is>
struct expander<T, std::index_sequence<Is...>> {
  template <typename E, std::size_t>
  using elem = E;

  using type = std::tuple<elem<T, Is>...>;
};

template <size_t N, typename T>
struct evaluated_tuple {
  using type = typename expander<T, std::make_index_sequence<N>>::type;
};

template <size_t order, typename T, size_t... Is>
inline typename evaluated_tuple<order, T>::type evaluate_into_tuple_impl(
    const state<T>& s, std::index_sequence<Is...>) {
  return std::make_tuple(evaluate<Is + 1>::get(s)...);
}

template <size_t order, typename T>
inline typename evaluated_tuple<order, T>::type evaluate_into_tuple(const state<T>& s) {
  return evaluate_into_tuple_impl<order>(s, std::make_index_sequence<order>{});
}

}  // namespace detail

template <typename A, typename B>
inline detail::state<B> init(const A& eccentricity, const B& mean_anomaly,
                             const B& eccentric_anomaly) {
  auto sincos = math::sincos(eccentric_anomaly);
  auto ecc_sin = B(eccentricity) * sincos.first;
  auto f0 = eccentric_anomaly - ecc_sin - mean_anomaly;
  return {f0, ecc_sin, B(eccentricity) * sincos.second};
}

template <size_t order, typename T>
inline T step(const detail::state<T>& state) {
  auto args = detail::evaluate_into_tuple<order>(state);
  return detail::householder_impl(state.f0, args, std::make_index_sequence<order>{});
}

}  // namespace householder
}  // namespace kepler

#endif