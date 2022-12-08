#ifndef KEPLER_MATH_HPP
#define KEPLER_MATH_HPP

#include <tuple>
#include <utility>

#include "simd.hpp"
#include "utils.hpp"

namespace kepler {
namespace math {

template <typename T>
inline T fma(const T& a, const T& b, const T& c) {
  return a * b + c;
}

template <typename A, typename T>
inline xs::batch<T, A> fma(const xs::batch<T, A>& a, const xs::batch<T, A>& b,
                           const xs::batch<T, A>& c) {
  return xs::fma(a, b, c);
}

template <typename T>
inline T fnma(const T& a, const T& b, const T& c) {
  return c - a * b;
}

template <typename A, typename T>
inline xs::batch<T, A> fnma(const xs::batch<T, A>& a, const xs::batch<T, A>& b,
                            const xs::batch<T, A>& c) {
  return xs::fnma(a, b, c);
}

template <typename T>
inline T horner_dynamic(const T&, const T& c1) {
  return c1;
}

template <typename T, typename... Args>
inline T horner_dynamic(const T& x, const T& c1, Args... c2) {
  return fma<T>(x, horner_dynamic(x, c2...), c1);
}

namespace detail {

template <std::uint64_t c>
struct coef {
  template <typename T>
  static inline T value() {
    return T(coef<c>::value<typename T::value_type>());
  }

  template <>
  static inline float value<float>() {
    return bit_cast<float>((uint32_t)c);
  }

  template <>
  static inline double value<double>() {
    return bit_cast<double>((uint64_t)c);
  }
};

}  // namespace detail

template <typename T, std::uint64_t c1>
inline T horner_static(const T&) noexcept {
  return detail::coef<c1>::template value<T>();
}

template <class T, std::uint64_t c1, std::uint64_t c2, std::uint64_t... args>
inline T horner_static(const T& x) noexcept {
  return fma<T>(x, horner_static<T, c2, args...>(x), detail::coef<c1>::template value<T>());
}

namespace detail {

template <typename T>
inline T short_sin_eval(const T&);

template <>
inline float short_sin_eval<float>(const float& x) {
  auto x2 = x * x;
  return x * horner_static<float, 0x3f800000, 0x3e2aaaab, 0x3c088889, 0x39500d01, 0x3638ef1d,
                           0x32d7322b, 0x2f309231, 0x2b573f9f>(-x2);
}
template <>
inline double short_sin_eval<double>(const double& x) {
  auto x2 = x * x;
  return x * horner_static<double, 0x3ff0000000000000, 0x3fc5555555555555, 0x3f81111111111111,
                           0x3f2a01a01a01a01a, 0x3ec71de3a556c734, 0x3e5ae64567f544e4,
                           0x3de6124613a86d09, 0x3d6ae7f3e733b81f>(-x2);
}

template <typename A>
inline xs::batch<float, A> short_sin_eval(const xs::batch<float, A>& x) {
  using B = xs::batch<float, A>;
  auto x2 = x * x;
  return x * horner_static<B, 0x3f800000, 0x3e2aaaab, 0x3c088889, 0x39500d01, 0x3638ef1d,
                           0x32d7322b, 0x2f309231, 0x2b573f9f>(-x2);
}

template <typename A>
inline xs::batch<double, A> short_sin_eval(const xs::batch<double, A>& x) {
  using B = xs::batch<double, A>;
  auto x2 = x * x;
  return x * horner_static<B, 0x3ff0000000000000, 0x3fc5555555555555, 0x3f81111111111111,
                           0x3f2a01a01a01a01a, 0x3ec71de3a556c734, 0x3e5ae64567f544e4,
                           0x3de6124613a86d09, 0x3d6ae7f3e733b81f>(-x2);
}

}  // namespace detail

// NOTE: only valid for x in [0, pi]
template <typename T>
inline std::pair<T, T> sincos(const T& x) {
  T s, c;
  if (x < constants::pio4<T>()) {
    s = detail::short_sin_eval(x);
    c = std::sqrt(1 - s * s);
  } else if (x > constants::threepio4<T>()) {
    s = detail::short_sin_eval(constants::pi<T>() - x);
    c = -std::sqrt(1 - s * s);
  } else {
    c = detail::short_sin_eval(constants::pio2<T>() - x);
    s = std::sqrt(1 - c * c);
  }
  return std::make_pair(s, c);
}

template <typename A, typename T>
inline std::pair<xs::batch<T, A>, xs::batch<T, A>> sincos(const xs::batch<T, A>& x) {
  return xs::sincos(x);
}

}  // namespace math
}  // namespace kepler

#endif
