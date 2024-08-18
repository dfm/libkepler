#ifndef KEPLER_MATH_HPP
#define KEPLER_MATH_HPP

#include <tuple>
#include <utility>

#include "kepler/kepler/utils.hpp"
#include "xsimd/xsimd.hpp"

namespace kepler {
namespace math {

namespace xs = xsimd;

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

template <typename T>
struct value_type {
  using type = typename T::value_type;
};

template <>
struct value_type<float> {
  using type = float;
};

template <>
struct value_type<double> {
  using type = double;
};

template <typename T>
struct as_unsigned_integer_t {
  using type = typename as_unsigned_integer_t<value_type<T>>::type;
};

template <>
struct as_unsigned_integer_t<float> {
  using type = std::uint32_t;
};

template <>
struct as_unsigned_integer_t<double> {
  using type = std::uint64_t;
};

template <typename T, std::uint64_t c>
inline T coef() {
  using V = typename value_type<T>::type;
  return T(bit_cast<V>((typename as_unsigned_integer_t<V>::type)(c)));
}

}  // namespace detail

template <typename T, std::uint64_t c1>
inline T horner_static(const T&) noexcept {
  return detail::coef<T, c1>();
}

template <class T, std::uint64_t c1, std::uint64_t c2, std::uint64_t... args>
inline T horner_static(const T& x) noexcept {
  return fma<T>(x, horner_static<T, c2, args...>(x), detail::coef<T, c1>());
}

namespace detail {

template <typename T>
static inline T cos_eval(const T&);
template <typename T>
static inline T sin_eval(const T&, const T&);

template <>
inline float cos_eval<float>(const float& z) {
  float y = horner_static<float, 0x3d2aaaa5, 0xbab60619, 0x37ccf5ce>(z);
  return 1. + fma(z, -0.5f, y * z * z);
}

template <>
inline float sin_eval<float>(const float& z, const float& x) {
  float y = horner_static<float, 0xbe2aaaa2, 0x3c08839d, 0xb94ca1f9>(z);
  return fma(y * z, x, x);
}

template <>
inline double cos_eval<double>(const double& z) {
  double y = horner_static<double, 0x3fe0000000000000ull, 0xbfa5555555555551ull,
                           0x3f56c16c16c15d47ull, 0xbefa01a019ddbcd9ull, 0x3e927e4f8e06d9a5ull,
                           0xbe21eea7c1e514d4ull, 0x3da8ff831ad9b219ull>(z);
  return 1. - y * z;
}

template <>
inline double sin_eval<double>(const double& z, const double& x) {
  double y =
      horner_static<double, 0xbfc5555555555548ull, 0x3f8111111110f7d0ull, 0xbf2a01a019bfdf03ull,
                    0x3ec71de3567d4896ull, 0xbe5ae5e5a9291691ull, 0x3de5d8fd1fcf0ec1ull>(z);
  return fma(y * z, x, x);
}

template <typename T>
inline unsigned short_reduce(const T& x, T& xr) {
  if (x < constants::pio4<T>()) {
    xr = x;
    return 0;
  } else if (x < constants::threepio4<T>()) {
    xr = x - constants::pio2<T>();
    return 1;
  }
  xr = x - constants::pi<T>();
  return 2;
}

}  // namespace detail

// NOTE: only valid for x in [0, pi], with a tiny bit of leeway...
template <typename T>
inline std::pair<T, T> sincos(const T& x) {
  const auto sgn = std::signbit(x) ? T(-1.) : T(1.);
  auto xr = std::abs(x);
  const auto n = detail::short_reduce(xr, xr);
  const auto z = xr * xr;
  const auto se = detail::sin_eval(z, xr);
  const auto ce = detail::cos_eval(z);
  if (n == 0) {
    return std::make_pair(sgn * se, ce);
  } else if (n == 1) {
    return std::make_pair(sgn * ce, -se);
  }
  return std::make_pair(sgn * -se, -ce);
}

template <typename A, typename T>
inline std::pair<xs::batch<T, A>, xs::batch<T, A>> sincos(const xs::batch<T, A>& x) {
  return xs::sincos(x);
}

}  // namespace math
}  // namespace kepler

#endif
