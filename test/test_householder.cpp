#include <cmath>

#include "./test_utils.hpp"

template <typename T>
inline T newton_reference(T e, T M, T E) {
  T fi = E - e * std::sin(E) - M;
  T fip = 1 - e * std::cos(E);
  return -fi / fip;
}

TEMPLATE_TEST_CASE("Newton", "[householder]", double, float) {
  using T = TestType;
  const T abs_tol = default_abs<TestType>::value;
  const std::size_t size = 100;
  const T e = 0.65;
  const T M = 0.123;
  for (std::size_t n = 0; n < size; ++n) {
    auto E = kepler::constants::pi<T>() * n / T(size - 1);
    auto expect = newton_reference(e, M, E);
    auto state = kepler::householder::householder<1>::init(e, M, E);
    auto calc = kepler::householder::householder<1>::step(state, e);
    REQUIRE_THAT(calc, WithinAbs(expect, abs_tol));
  }
}

template <typename T>
inline T halley_reference(T e, T M, T E) {
  T fi = E - e * std::sin(E) - M;
  T fip = 1 - e * std::cos(E);
  T fipp = e * std::sin(E);
  T fippp = 1 - fip;
  T d1 = -fi / fip;
  d1 = -fi / (fip + d1 * fipp / 2.0);
  d1 = -fi / (fip + d1 * fipp / 2.0 + d1 * d1 * fippp / 6.0);
  return d1;
}

TEMPLATE_TEST_CASE("Halley", "[householder]", double, float) {
  using T = TestType;
  const T abs_tol = default_abs<TestType>::value;
  const std::size_t size = 100;
  const T e = 0.65;
  const T M = 0.123;
  for (std::size_t n = 0; n < size; ++n) {
    auto E = kepler::constants::pi<T>() * n / T(size - 1);
    auto expect = halley_reference(e, M, E);
    auto state = kepler::householder::householder<3>::init(e, M, E);
    auto calc = kepler::householder::householder<3>::step(state, e);
    REQUIRE_THAT(calc, WithinAbs(expect, abs_tol));
  }
}
