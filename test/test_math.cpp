#include <cmath>

#include "./test_utils.hpp"
#include "kepler/kepler/constants.hpp"
#include "kepler/kepler/math.hpp"

using namespace kepler;

TEMPLATE_TEST_CASE("Short sine", "[math]", double, float) {
  using T = TestType;
  const T abs_tol = default_abs<TestType>::value;
  const std::size_t size = 10000;
  for (std::size_t n = 0; n < size; ++n) {
    auto x = constants::pi<T>() * n / T(size - 1);
    auto calc = math::sincos(x);
    REQUIRE_THAT(calc.first, WithinAbs(std::sin(x), abs_tol));
    REQUIRE_THAT(calc.second, WithinAbs(std::cos(x), abs_tol));
  }
}

TEMPLATE_TEST_CASE("Short sine (SIMD)", "[math][simd]", double, float) {
  using T = TestType;
  const T abs_tol = default_abs<TestType>::value;
  const std::size_t size = 10000;

  for (std::size_t n = 0; n < size; ++n) {
    auto x = constants::pi<T>() * n / T(size - 1);
    auto calc = math::sincos(xs::batch<T>(x));
    REQUIRE_THAT(calc.first.get(0), WithinAbs(std::sin(x), abs_tol));
    REQUIRE_THAT(calc.second.get(0), WithinAbs(std::cos(x), abs_tol));
  }
}
