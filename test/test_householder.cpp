#include <cmath>

#include "./test_utils.hpp"

TEMPLATE_TEST_CASE("Evaluate", "[householder]", double, float) {
  using namespace kepler::householder::detail;
  using T = TestType;
  const T abs_tol = default_abs<TestType>::value;
  const T e = 0.65;
  const T M = 0.123;
  const T E = 0.456;
  auto state = kepler::householder::init(e, M, E);

  REQUIRE_THAT(evaluate<1>::get(state), WithinAbs(1 - e * std::cos(E), abs_tol));
  REQUIRE_THAT(evaluate<2>::get(state), WithinAbs(e * std::sin(E), abs_tol));
  REQUIRE_THAT(evaluate<3>::get(state), WithinAbs(e * std::cos(E), abs_tol));
  REQUIRE_THAT(evaluate<4>::get(state), WithinAbs(-e * std::sin(E), abs_tol));
  REQUIRE_THAT(evaluate<5>::get(state), WithinAbs(-e * std::cos(E), abs_tol));
  REQUIRE_THAT(evaluate<6>::get(state), WithinAbs(e * std::sin(E), abs_tol));
  REQUIRE_THAT(evaluate<7>::get(state), WithinAbs(e * std::cos(E), abs_tol));
  REQUIRE_THAT(evaluate<8>::get(state), WithinAbs(-e * std::sin(E), abs_tol));
  REQUIRE_THAT(evaluate<9>::get(state), WithinAbs(-e * std::cos(E), abs_tol));
}

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
    auto state = kepler::householder::init(e, M, E);
    auto calc = kepler::householder::step<1>(state);
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
    auto state = kepler::householder::init(e, M, E);
    auto calc = kepler::householder::step<3>(state);
    REQUIRE_THAT(calc, WithinAbs(expect, abs_tol));
  }
}

template <typename T>
inline T fourth_reference(T e, T M, T E) {
  T fi = E - e * std::sin(E) - M;
  T fip = 1 - e * std::cos(E);
  T fipp = e * std::sin(E);
  T fippp = 1 - fip;
  T fipppp = -fipp;
  T d1 = -fi / fip;
  d1 = -fi / (fip + d1 * fipp / 2.0);
  d1 = -fi / (fip + d1 * fipp / 2.0 + d1 * d1 * fippp / 6.0);
  d1 = -fi / (fip + d1 * fipp / 2.0 + d1 * d1 * fippp / 6.0 + d1 * d1 * d1 * fipppp / 24.0);
  return d1;
}

TEMPLATE_TEST_CASE("Fourth", "[householder]", double, float) {
  using T = TestType;
  const T abs_tol = default_abs<TestType>::value;
  const std::size_t size = 100;
  const T e = 0.65;
  const T M = 0.123;
  for (std::size_t n = 0; n < size; ++n) {
    auto E = kepler::constants::pi<T>() * n / T(size - 1);
    auto expect = fourth_reference(e, M, E);
    auto state = kepler::householder::init(e, M, E);
    auto calc = kepler::householder::step<4>(state);
    REQUIRE_THAT(calc, WithinAbs(expect, abs_tol));
  }
}
