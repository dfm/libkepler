#ifndef KEPLER_MATH_HPP
#define KEPLER_MATH_HPP

#include <tuple>
#include <utility>

#include "./simd.hpp"

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
inline T horner(const T&, const T& c1) {
  return c1;
}

template <typename T, typename... Args>
inline T horner(const T& x, const T& c1, Args... c2) {
  return fma<T>(horner(x, c2...), x, c1);
}

namespace detail {

template <typename T>
inline T short_sin_eval(const T& x) {
  auto x2 = x * x;
  return x * horner(-x2, T(1), constants::shortsin1<T>(), constants::shortsin2<T>(),
                    constants::shortsin3<T>(), constants::shortsin4<T>(),
                    constants::shortsin5<T>(), constants::shortsin6<T>(),
                    constants::shortsin7<T>());
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
