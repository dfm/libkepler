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
