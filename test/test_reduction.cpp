#include "./test_utils.hpp"

TEST_CASE("Range reduction", "[reduction]") {
  const double abs_tol = 5e-15;
  double xr;

  REQUIRE(kepler::range_reduce(0.0, xr) == false);
  REQUIRE_THAT(xr, WithinAbs(0.0, abs_tol));

  REQUIRE(kepler::range_reduce(kepler::constants::pi<double>(), xr) == true);
  REQUIRE_THAT(xr, WithinAbs(kepler::constants::pi<double>(), 1e-15));

  kepler::range_reduce(kepler::constants::twopi<double>(), xr);
  REQUIRE_THAT(xr, WithinAbs(0.0, abs_tol));

  REQUIRE(kepler::range_reduce(100 * kepler::constants::twopi<double>() - 1e-8, xr) == true);
  REQUIRE_THAT(xr, WithinAbs(1e-8, abs_tol));
}

TEST_CASE("Range reduction (SIMD)", "[reduction][simd]") {
  typedef double T;
  using B = xs::batch<T>;
  const T abs_tol = 5e-15;
  B xr;

  REQUIRE(!xs::any(kepler::range_reduce(B(0.0), xr)));
  REQUIRE_THAT(xr.get(0), WithinAbs(0.0, abs_tol));

  REQUIRE(xs::all(kepler::range_reduce(kepler::constants::pi<B>(), xr)));
  REQUIRE_THAT(xr.get(0), WithinAbs(kepler::constants::pi<T>(), 1e-15));

  kepler::range_reduce(kepler::constants::twopi<B>(), xr);
  REQUIRE_THAT(xr.get(0), WithinAbs(0.0, abs_tol));

  REQUIRE(xs::all(kepler::range_reduce(B(100.) * kepler::constants::twopi<B>() - B(1e-8), xr)));
  REQUIRE_THAT(xr.get(0), WithinAbs(1e-8, abs_tol));
}

// https://stackoverflow.com/questions/42792939/implementation-of-sinpi-and-cospi-using-standard-c-math-library/42792940#42792940
template <typename T>
struct int_type {};

template <>
struct int_type<float> {
  using type = int32_t;
};

template <>
struct int_type<double> {
  using type = int64_t;
};

template <typename T>
inline typename int_type<T>::type reduce(T& x) noexcept;

template <>
inline int32_t reduce<float>(float& x) noexcept {
  auto xz = x * 0.0f;
  x = (fabsf(x) < 0x1.0p24f) ? x : xz;
  auto r = nearbyintf(x + x);
  x = fmaf(-0.5f, r, x);
  return int32_t(r);
}

template <>
inline int64_t reduce<double>(double& x) noexcept {
  auto xz = x * 0.0;
  x = (fabs(x) < 9.0071992547409920e+15) ? x : xz;
  auto r = nearbyint(x + x);
  x = fma(-0.5f, r, x);
  return int64_t(r);
}

TEST_CASE("Faster range reduction") {
  double x;
  int64_t q;

  x = 1000. + 0.1;
  q = reduce(x);
  REQUIRE_THAT(x, WithinAbs(1e-12, 0.1));

  x = 1000. + 0.6;
  q = reduce(x);
  REQUIRE(q & 1);
  REQUIRE_THAT(x, WithinAbs(1e-12, 0.1));
}

// az = a * 0.0; // must be evaluated with IEEE-754 semantics
// /* for |a| >= 2**53, cospi(a) = 1.0, but cospi(Inf) = NaN */
// a = (fabs (a) < 9.0071992547409920e+15) ? a : az;  // 0x1.0p53
// /* reduce argument to primary approximation interval (-0.25, 0.25) */
// r = nearbyint (a + a); // must use IEEE-754 "to nearest" rounding
// i = (int64_t)r;
// t = fma (-0.5, r, a);
