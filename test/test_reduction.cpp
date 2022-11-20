#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "kepler/kepler.hpp"

using namespace Catch::Matchers;

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
