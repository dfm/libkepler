#ifndef _KEPLER_BENCHMARKS_SIMD_H_
#define _KEPLER_BENCHMARKS_SIMD_H_

#include "xsimd/xsimd.hpp"

namespace xs = xsimd;

namespace kepler_benchmarks {
namespace simd {

#define KB_SIMD_DEFINE_CONSTANT(NAME, SINGLE, DOUBLE) \
  template <class T>                                  \
  inline T NAME() noexcept {                          \
    return T(NAME<typename T::value_type>());         \
  }                                                   \
  template <>                                         \
  inline float NAME<float>() noexcept {               \
    return xs::bit_cast<float>((uint32_t)SINGLE);     \
  }                                                   \
  template <>                                         \
  inline double NAME<double>() noexcept {             \
    return xs::bit_cast<double>((uint64_t)DOUBLE);    \
  }

KB_SIMD_DEFINE_CONSTANT(twopi, 0x40c90fdb, 0x401921fb54442d18)
KB_SIMD_DEFINE_CONSTANT(sixth, 0x3e2aaaab, 0x3fc5555555555555)
KB_SIMD_DEFINE_CONSTANT(twentieth, 0x3d2aaaab, 0x3fa5555555555555)
KB_SIMD_DEFINE_CONSTANT(factor1, 0x40f4da39, 0x401e9b471164c596)
KB_SIMD_DEFINE_CONSTANT(factor2, 0x3fa6450f, 0x3ff4c8a1d518acbd)

#undef KB_SIMD_DEFINE_CONSTANT

template <class B, uint64_t c0>
inline B poly(const B& x) noexcept {
  return B(1.) - xs::kernel::detail::coef<B, c0>() * x;
}

template <class B, uint64_t c0, uint64_t c1, uint64_t... args>
inline B poly(const B& x) noexcept {
  return xs::fnma(xs::kernel::detail::coef<B, c0>() * x, poly<B, c1, args...>(x), B(1.));
}

// sin_reduc(x, x^2) = x - sin(x); for x in [0, pi/4)
template <class A, class T>
inline xs::batch<T, A> sin_reduc_eval(const xs::batch<T, A>& x,
                                      const xs::batch<T, A>& x2) noexcept;

template <class A>
inline xs::batch<float, A> sin_reduc_eval(const xs::batch<float, A>& x,
                                          const xs::batch<float, A>& x2) noexcept {
  auto s = poly<xs::batch<float, A>, 0x3d4ccccd, 0x3cc30c31, 0x3c638e39, 0x3c14f209, 0x3bd20d21,
                0x3b9c09c1, 0x3b70f0f1, 0x3b3fa030, 0x3b1c09c1>(x2);
  s *= x * x2 * xs::kernel::detail::coef<xs::batch<float, A>, 0x3e2aaaab>();
  return s;
}

template <class A>
inline xs::batch<double, A> sin_reduc_eval(const xs::batch<double, A>& x,
                                           const xs::batch<double, A>& x2) noexcept {
  auto s = poly<xs::batch<double, A>, 0x3fa999999999999a, 0x3f98618618618618, 0x3f8c71c71c71c71c,
                0x3f829e4129e4129e, 0x3f7a41a41a41a41a, 0x3f73813813813814, 0x3f6e1e1e1e1e1e1e,
                0x3f67f405fd017f40, 0x3f63813813813814>(x2);
  s *= x * x2 * xs::kernel::detail::coef<xs::batch<double, A>, 0x3fc5555555555555>();
  return s;
}

// cos_reduc(x^2) = 1 - cos(x); for x in [0, pi/4)
template <class A, class T>
inline xs::batch<T, A> cos_reduc_eval(const xs::batch<T, A>& x2) noexcept;

template <class A>
inline xs::batch<float, A> cos_reduc_eval(const xs::batch<float, A>& x2) noexcept {
  auto c = poly<xs::batch<float, A>, 0x3daaaaab, 0x3d088889, 0x3c924925, 0x3c360b61, 0x3bf83e10,
                0x3bb40b41, 0x3b888889, 0x3b562b81, 0x3b2c7692>(x2);
  c *= x2 * xs::kernel::detail::coef<xs::batch<float, A>, 0x3f000000>();
  return c;
}

template <class A>
inline xs::batch<double, A> cos_reduc_eval(const xs::batch<double, A>& x2) noexcept {
  auto c = poly<xs::batch<double, A>, 0x3fb5555555555555, 0x3fa1111111111111, 0x3f92492492492492,
                0x3f86c16c16c16c17, 0x3f7f07c1f07c1f08, 0x3f76816816816817, 0x3f71111111111111,
                0x3f6ac5701ac5701b, 0x3f658ed2308158ed>(x2);
  c *= x2 * xs::kernel::detail::coef<xs::batch<double, A>, 0x3fe0000000000000>();
  return c;
}

// x - sin(x) and 1 - cos(x) for x in [0, pi)
template <class B>
inline void sin_cos_reduc(const B& x, B& sr, B& cr) noexcept {
  auto big = x >= xs::constants::pio2<B>();
  auto u = xs::select(big, xs::constants::pi<B>() - x, x);

  auto mid = u >= xs::constants::pio2<B>();
  auto v = xs::select(mid, xs::constants::pi<B>() - u, u);
  auto v2 = v * v;

  sr = sin_reduc_eval(v, v2);
  cr = cos_reduc_eval(v2);

  sr = xs::select(mid, u + cr - B(1.), sr);
  cr = xs::select(mid, u + sr + B(1. - xs::constants::pio2<typename B::value_type>()), cr);

  sr = xs::select(big, xs::fma(B(2.), x, sr) - xs::constants::pi<B>(), sr);
  cr = xs::select(big, B(2.) - cr, cr);
}

// Numerically stable range reduction of x into [0, pi). Returns true for x in [pi, 2*pi)
template <class A, class T>
inline xs::batch_bool<T, A> reduce(xs::batch<T, A> const& x, xs::batch<T, A>& xr) noexcept {
  using B = xs::batch<T, A>;
  auto quad = xs::kernel::detail::trigo_reducer<B>::reduce(x, xr);
  xr = xs::fma(quad, xs::constants::pio2<B>(), xr);
  auto lo = xr < B(0.);
  auto hi = xr >= xs::constants::pi<B>();
  xr = xs::select(hi, twopi<B>() - xr, xs::abs(xr));
  return hi || lo;
}

template <class A, class T>
inline void starter(const T& ecc, const T& ome, const xs::batch<T, A>& mean_anom,
                    xs::batch<T, A>& ecc_anom) noexcept {
  using B = xs::batch<T, A>;

  auto Bome = B(ome);
  auto m2 = mean_anom * mean_anom;

  // alpha = FACTOR1 + FACTOR2 * (M_PI - M) / (1 + ecc)
  auto alpha =
      xs::fma(B(factor2<T>() / (1. + ecc)), xs::constants::pi<B>() - mean_anom, factor1<B>());

  // d = 3 * ome + alpha * ecc
  auto d = xs::fma(B(ecc), alpha, B(3. * ome));
  alpha *= d;

  // r = (3 * alpha * (d - ome) + M2) * M
  auto r = mean_anom * xs::fma(B(3.) * alpha, d - Bome, m2);

  // q = 2 * alphad * ome - M2
  auto q = xs::fms(B(2.) * alpha, Bome, m2);
  auto q2 = q * q;

  // w = cbrt(std::abs(r) + sqrt(q2 * q + r * r))
  auto w = xs::cbrt(xs::abs(r) + xs::sqrt(xs::fma(q2, q, r * r)));
  w *= w;

  // E = (2 * r * w / (w * w + w * q + q2) + M) / d
  auto denom = xs::fma(w, w + q, q2);
  ecc_anom = xs::fma(B(2.) * r / denom, w, mean_anom) / d;
}

template <class B>
inline B newton_step(const B& f0, const B& f1, const B& f2, const B& f3) noexcept {
  // d_3 = -f_0 / (f_1 - 0.5 * f_0 * f_2 / f_1)
  auto d3 = f0 / xs::fms(B(0.5) * f0 / f1, f2, f1);

  // d_4 = -f_0 / (f_1 + 0.5 * d_3 * f_2 + (d_3 * d_3) * f_3 / 6)
  auto d4 = f0 / xs::fnms(d3, xs::fma(sixth<B>() * d3, f3, B(0.5) * f2), f1);

  // dE = -f_0 / (f_1 + 0.5 * d_4 * f_2 + d_4 * d_4 * f_3 / 6 - d_42 * d_4 * f_2 / 24)
  return f0 /
         xs::fnms(d4, xs::fma(d4, xs::fnma(twentieth<B>() * d4, f2, sixth<B>() * f3), B(0.5) * f2),
                  f1);
}

template <class A, class T>
inline void refine(const T& ecc, const T& ome, const xs::batch<T, A>& mean_anom,
                   xs::batch<T, A>& ecc_anom) noexcept {
  using B = xs::batch<T, A>;
  auto Bome = B(ome);
  auto Becc = B(ecc);
  B sr, cr;
  sin_cos_reduc(ecc_anom, sr, cr);
  sr *= Becc;
  auto f0 = sr + xs::fms(ecc_anom, Bome, mean_anom);
  auto f1 = xs::fma(Becc, cr, Bome);
  auto f2 = xs::fms(Becc, ecc_anom, sr);
  auto f3 = B(1.) - f1;
  ecc_anom += newton_step(f0, f1, f2, f3);
}

template <class A, class T>
inline xs::batch<T, A> solve_one_batch(const T& ecc, const xs::batch<T, A>& mean_anom) noexcept {
  using B = xs::batch<T, A>;
  auto ome = 1. - ecc;
  B mean_anom_reduc, ecc_anom;
  auto flag = reduce(mean_anom, mean_anom_reduc);
  starter(ecc, ome, mean_anom_reduc, ecc_anom);
  refine(ecc, ome, mean_anom_reduc, ecc_anom);
  return xs::select(flag, twopi<B>() - ecc_anom, ecc_anom);
}

template <class T, class Tag = xs::unaligned_mode>
inline void solve(std::size_t size, const T& ecc, const T* mean_anom, T* ecc_anom) noexcept {
  constexpr std::size_t simd_size = xs::simd_type<T>::size;
  std::size_t vec_size = size - size % simd_size;

  for (std::size_t i = 0; i < vec_size; i += simd_size) {
    auto m = xs::load(&mean_anom[i], Tag());
    auto result = solve_one_batch(ecc, m);
    xs::store(&ecc_anom[i], result, Tag());
  }

  if (vec_size < size) {
    alignas(xs::batch<T>::arch_type::alignment()) std::array<T, simd_size> buffer;
    buffer.fill(T(0.));
    for (std::size_t i = vec_size + 1; i < size; ++i) {
      buffer[i - vec_size] = mean_anom[i];
    }
    auto m = xs::load_aligned(&buffer[0]);
    auto result = solve_one_batch(ecc, m);
    result.store_aligned(&buffer[vec_size]);
    for (std::size_t i = vec_size + 1; i < size; ++i) {
      ecc_anom[i] = buffer[i - vec_size];
    }
  }
}

}  // namespace simd
}  // namespace kepler_benchmarks

#endif
