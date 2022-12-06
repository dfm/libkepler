#ifndef KEPLER_REDUCTION_HPP
#define KEPLER_REDUCTION_HPP

#include <cmath>
#include <cstdint>
#include <limits>

#include "constants.hpp"
#include "simd.hpp"

namespace kepler {

namespace detail {

// The following functions perform argument range reductions for all the solvers
// implemented in this library. The implementation is based closely on the
// trigonometric reduction algorithms used by the BSD-licensed xsimd library.
//
// Source: xsimd/arch/generic/xsimd_generic_trigo.hpp
//
// The xsimd license is:
//
/****************************************************************************
 * Copyright (c) Johan Mabille, Sylvain Corlay, Wolf Vollprecht and         *
 * Martin Renou                                                             *
 * Copyright (c) QuantStack                                                 *
 * Copyright (c) Serge Guelton                                              *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 ****************************************************************************/

template <typename T>
inline T quadrant(const T& x) noexcept {
  return x & T(3);
}

template <>
inline float quadrant(const float& x) noexcept {
  return static_cast<float>(quadrant(static_cast<int>(x)));
}

template <>
inline double quadrant(const double& x) noexcept {
  double a = x * 0.25;
  return 4. * (a - floor(a));
}

template <typename T>
inline T trig_reduce(const T& x, T& xr) noexcept {
  if (x <= constants::pio4<T>()) {
    xr = x;
    return T(0.);
  } else if (x <= constants::pio2<T>()) {
    xr = x - constants::pio2_1<T>();
    xr -= constants::pio2_2<T>();
    xr -= constants::pio2_3<T>();
    return T(1.);
  } else if (x <= constants::twentypi<T>()) {
    T xi = nearbyint(x * constants::twoopi<T>());
    xr = x - xi * constants::pio2_1<T>();
    xr -= xi * constants::pio2_2<T>();
    xr -= xi * constants::pio2_3<T>();
    return quadrant(xi);
  } else if (x <= constants::mediumpi<T>()) {
    T fn = nearbyint(x * constants::twoopi<T>());
    T r = x - fn * constants::pio2_1<T>();
    T w = fn * constants::pio2_1t<T>();
    T t = r;
    w = fn * constants::pio2_2<T>();
    r = t - w;
    w = fn * constants::pio2_2t<T>() - ((t - r) - w);
    t = r;
    w = fn * constants::pio2_3<T>();
    r = t - w;
    w = fn * constants::pio2_3t<T>() - ((t - r) - w);
    xr = r - w;
    return quadrant(fn);
  } else {
    T n;
    if (x == std::numeric_limits<T>::infinity()) {
      n = T(0.);
      xr = std::numeric_limits<T>::quiet_NaN();
    } else {
      double y[2];
      std::int32_t n_ = xsimd::detail::__ieee754_rem_pio2(x, y);
      n = T(n_ & 3);
      xr = T(y[0]);
    }
    return n;
  }
}

}  // namespace detail

template <typename T>
inline bool range_reduce(const T& x, T& xr) noexcept {
  auto quad = detail::trig_reduce(x, xr);
  xr += quad * constants::pio2<T>();
  if (xr < T(0.)) {
    xr = -xr;
    return true;
  } else if (xr >= constants::pi<T>()) {
    xr = constants::twopi<T>() - xr;
    return true;
  }
  return false;
}

template <typename A, typename T>
inline xs::batch_bool<T, A> range_reduce(xs::batch<T, A> const& x, xs::batch<T, A>& xr) noexcept {
  using B = xs::batch<T, A>;
  auto quad = xs::kernel::detail::trigo_reducer<B>::reduce(x, xr);
  xr = xs::fma(quad, constants::pio2<B>(), xr);
  auto lo = xr < B(0.);
  auto hi = xr >= constants::pi<B>();
  xr = xs::select(hi, constants::twopi<B>() - xr, xs::abs(xr));
  return hi | lo;
}

}  // namespace kepler
#endif
