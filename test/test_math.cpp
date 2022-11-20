#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "kepler/kepler.hpp"

using namespace Catch::Matchers;

TEST_CASE("Sin and cos reduced", "[math]") {
  using T = double;
  const T abs_tol = 1e-15;
  const std::size_t size = 1000;

  for (size_t n = 0; n < size; ++n) {
    auto x = kepler::constants::pi<T>() * n / T(size - 1);
    T sr, cr;
    kepler::detail::sin_cos_reduc(x, sr, cr);
    REQUIRE_THAT(sr, WithinAbs(x - std::sin(x), abs_tol));
    REQUIRE_THAT(cr, WithinAbs(1 - std::cos(x), abs_tol));
  }
}
