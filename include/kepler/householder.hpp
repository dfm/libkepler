///
/// This header provides the functions needed to construct root finding updates
/// using Householder's method to arbitrary order. This is a crucial element of
/// just about all Kepler solvers since most solvers require either repeated
/// iterations of Householder's method or a single iteration starting from a
/// high precision initial guess.
///
/// The is mainly used as follows:
///
/// ```cpp
/// auto state = kepler::householder::init(eccen, mean_anom, ecc_anom_guess);
/// auto dE = kepler::householder::step<ORDER>(state);
/// ```
///
/// where `ORDER` is any integer `>= 1`. To get a Newton step, use `ORDER = 1`,
/// and `ORDER = 3` gives Halley's method.
///
/// For iterative methods, the convergence test will generally be executed using
/// `state.f0` which is the difference between `mean_anom` and the mean anomaly
/// computed from `eccen` and `ecc_anom_guess`.
///
/// The implementation details might seem a little convoluted, but the idea here
/// is to transparently support arbitrary values of `ORDER` without sacrificing
/// performance. With sensible compiler optimization settings, the methods here
/// produce machine code that is equivalent to the hand-written versions, thanks
/// to some gnarly template metaprogramming.
///

#ifndef KEPLER_HOUSEHOLDER_HPP
#define KEPLER_HOUSEHOLDER_HPP

#include <cstdint>

#include "./math.hpp"

namespace kepler {
namespace householder {
namespace detail {

/// The `state` contains the info needed to perform an update
template <typename T>
struct state {
  T f0;       ///< `E - e*sin(E) - M` at the current iteration
  T ecc_sin;  ///< `e * sin(E)` at the current iteration
  T ecc_cos;  ///< `e * cos(E)` at the current iteration
};

/// The following provides an interface for computing arbitrary order of
/// derivatives of `f0 = E - e*sin(E) - M` with respect to `E`. The key
/// realization is that for `n > 3`, the `n`th derivative of `f0` can be
/// computed simply as `+/- e * sin(E)` when `n` is even, and `+/- e * cos(E)`
/// otherwise. The sign in each case depends on the parity of `n/2`.
///
/// In practice, `evaluate<n>::get(state)` will return the `n`th derivative for
/// any `n >= 1`.
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

/// Next we have methods for evaluating factorials at compile time, since these
/// are used to normalize the terms in the Householder step. In practice,
/// `factorial<n>::value` evaluates to the compile time constant `n!`.
template <std::uint64_t n>
struct factorial {
  constexpr static std::uint64_t value = n * factorial<n - 1>::value;
};

template <>
struct factorial<0> {
  constexpr static std::uint64_t value = 1;
};

/// The `horner_packed` function provides an interface for evaluating a
/// polynomial using Horner's method with support for unpacking a tuple of
/// coefficients. This depends on the (very simple) `horner` function defined in
/// `math.hpp`.
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

/// Now we get into the real magic of this implementation. The `householder`
/// (and related) functions compute the Householder step for a given (arbitrary)
/// order. This should have been simple, but the problem is that the
/// implementation needs to be called recursively with all but the last
/// coefficient repeated at each lower order. This isn't possible to do with C++
/// parameter packs so instead we're packing and unpacking tuples of
/// coefficients. This seems like a bad idea, but all the standard compilers
/// seems to be clever enough to optimize out any overhead, producing sensible a
/// sensible implementation.
///
/// One key reference was:
/// https://aherrmann.github.io/programming/2016/02/28/unpacking-tuples-in-cpp14/
template <typename Tuple, size_t... Is>
constexpr auto tuple_init_impl(Tuple t, std::index_sequence<Is...>) {
  return std::make_tuple(std::get<Is>(t)...);
}

template <typename Tuple, typename Is = std::make_index_sequence<std::tuple_size<Tuple>{} - 1>>
constexpr auto tuple_init(Tuple t) {
  return tuple_init_impl(t, Is{});
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

/// The `evaluated_tuple` type is used for inferring the tuple type arguments
/// for a set of coefficients of a given order. This implementation is based on
/// the answer here: https://stackoverflow.com/a/66255219
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

/// The `evaluate_into_tuple` function is used to construct a tuple of
/// coefficients for a given order using the `evaluate` struct from above.
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