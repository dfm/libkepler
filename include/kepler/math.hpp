#ifndef KEPLER_MATH_HPP
#define KEPLER_MATH_HPP

#include "./constants.hpp"
#include "./utils.hpp"

namespace kepler {
namespace detail {

template <class T>
struct as_unsigned_integer : std::make_unsigned<T> {};

template <>
struct as_unsigned_integer<float> {
  using type = uint32_t;
};

template <>
struct as_unsigned_integer<double> {
  using type = uint64_t;
};

template <typename T>
using as_unsigned_integer_t = typename as_unsigned_integer<T>::type;

template <class T, uint64_t c>
inline T coef() noexcept {
  return bit_cast<T>(as_unsigned_integer_t<T>(c));
}

template <uint64_t c>
inline float coef() noexcept {
  return bit_cast<float>((uint32_t)c);
}

template <uint64_t c>
inline double coef() noexcept {
  return bit_cast<double>((uint64_t)c);
}

template <typename T, uint64_t c0>
inline T poly(const T& x) noexcept {
  return T(1.) - coef<T, c0>() * x;
}

template <typename T, uint64_t c0, uint64_t c1, uint64_t... args>
inline T poly(const T& x) noexcept {
  return T(1.) - coef<T, c0>() * x * poly<T, c1, args...>(x);
}

// sin_reduc(x, x^2) = x - sin(x); for x in [0, pi/4)
template <typename T>
inline T sin_reduc_eval(const T& x, const T& x2) noexcept;

template <>
inline float sin_reduc_eval<float>(const float& x, const float& x2) noexcept {
  auto s = poly<float, 0x3d4ccccd, 0x3cc30c31, 0x3c638e39, 0x3c14f209, 0x3bd20d21, 0x3b9c09c1,
                0x3b70f0f1, 0x3b3fa030, 0x3b1c09c1>(x2);
  s *= x * x2 * coef<float, 0x3e2aaaab>();
  return s;
}

template <>
inline double sin_reduc_eval<double>(const double& x, const double& x2) noexcept {
  auto s = poly<double, 0x3fa999999999999a, 0x3f98618618618618, 0x3f8c71c71c71c71c,
                0x3f829e4129e4129e, 0x3f7a41a41a41a41a, 0x3f73813813813814, 0x3f6e1e1e1e1e1e1e,
                0x3f67f405fd017f40, 0x3f63813813813814>(x2);
  s *= x * x2 * coef<double, 0x3fc5555555555555>();
  return s;
}

// cos_reduc(x^2) = 1 - cos(x); for x in [0, pi/4)
template <typename T>
inline T cos_reduc_eval(const T& x2) noexcept;

template <>
inline float cos_reduc_eval(const float& x2) noexcept {
  auto c = poly<float, 0x3daaaaab, 0x3d088889, 0x3c924925, 0x3c360b61, 0x3bf83e10, 0x3bb40b41,
                0x3b888889, 0x3b562b81, 0x3b2c7692>(x2);
  c *= x2 * coef<float, 0x3f000000>();
  return c;
}

template <>
inline double cos_reduc_eval(const double& x2) noexcept {
  auto c = poly<double, 0x3fb5555555555555, 0x3fa1111111111111, 0x3f92492492492492,
                0x3f86c16c16c16c17, 0x3f7f07c1f07c1f08, 0x3f76816816816817, 0x3f71111111111111,
                0x3f6ac5701ac5701b, 0x3f658ed2308158ed>(x2);
  c *= x2 * coef<double, 0x3fe0000000000000>();
  return c;
}

// x - sin(x) and 1 - cos(x) for x in [0, pi)
template <typename T>
inline void sin_cos_reduc(const T& x, T& sr, T& cr) noexcept {
  auto big = x >= constants::pio2<T>();
  auto u = big ? constants::pi<T>() - x : x;

  auto mid = u >= constants::pio2<T>();
  auto v = mid ? constants::pi<T>() - u : u;
  auto v2 = v * v;

  auto s = sin_reduc_eval(v, v2);
  auto c = cos_reduc_eval(v2);

  if (mid) {
    sr = u + c - T(1.);
    cr = u + s + T(1.) - constants::pio2<T>();
  } else {
    sr = s;
    cr = c;
  }
  if (big) {
    sr += T(2.) * x - constants::pi<T>();
    cr = T(2.) - cr;
  }
}

}  // namespace detail
}  // namespace kepler

#endif
