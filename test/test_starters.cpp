#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "kepler/constants.hpp"
#include "kepler/starters.hpp"

using namespace Catch::Matchers;

TEST_CASE("Mikkola starter", "[starters]") {
  using T = double;
  const T abs_tol = 1e-13;
  const T rel_tol = 0.002;  // This is the expected relative tolerance from the paper.
  const size_t ecc_size = 10;
  const size_t anom_size = 100;
  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    const kepler::starters::mikkola starter(eccentricity);
    for (size_t m = 0; m < anom_size; ++m) {
      const T ecc_anom_expect = kepler::constants::pi<T>() * m / T(anom_size - 1);
      auto mean_anomaly = ecc_anom_expect - eccentricity * std::sin(ecc_anom_expect);
      auto ecc_anom = starter.start(mean_anomaly);
      if (m == 0) {
        REQUIRE_THAT(ecc_anom, WithinAbs(ecc_anom_expect, abs_tol));
      } else {
        REQUIRE_THAT(ecc_anom, WithinRel(ecc_anom_expect, rel_tol));
      }
    }
  }
}

TEST_CASE("Markley starter", "[starters]") {
  using T = double;
  const T abs_tol = 1e-13;
  const T rel_tol = 5e-4;
  const size_t ecc_size = 10;
  const size_t anom_size = 100;
  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    const kepler::starters::markley starter(eccentricity);
    for (size_t m = 0; m < anom_size; ++m) {
      const T ecc_anom_expect = kepler::constants::pi<T>() * m / T(anom_size - 1);
      auto mean_anomaly = ecc_anom_expect - eccentricity * std::sin(ecc_anom_expect);
      auto ecc_anom = starter.start(mean_anomaly);
      if (m == 0) {
        REQUIRE_THAT(ecc_anom, WithinAbs(ecc_anom_expect, abs_tol));
      } else {
        REQUIRE_THAT(ecc_anom, WithinRel(ecc_anom_expect, rel_tol));
      }
    }
  }
}

TEST_CASE("RPP17/B21 singular corner", "[starters]") {
  using T = double;
  const T abs_tol = 1e-13;
  const T rel_tol = 1e-6;
  const size_t ecc_size = 10;
  const size_t anom_size = 50;
  const T ecc_min = 0.8;
  const T ecc_anom_max = 0.1;
  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = ecc_min + (1. - ecc_min) * n / T(ecc_size);
    const kepler::starters::rppb starter(eccentricity);
    for (size_t m = 0; m < anom_size; ++m) {
      const T ecc_anom_expect = ecc_anom_max * m / T(anom_size - 1);
      auto mean_anomaly = ecc_anom_expect - eccentricity * std::sin(ecc_anom_expect);
      auto ecc_anom = starter.singular(mean_anomaly);
      if (m == 0) {
        REQUIRE_THAT(ecc_anom, WithinAbs(ecc_anom_expect, abs_tol));
      } else {
        REQUIRE_THAT(ecc_anom, WithinRel(ecc_anom_expect, rel_tol));
      }
    }
  }
}

TEST_CASE("RPP17/B21 starter", "[starters]") {
  using T = double;
  const T abs_tol = 1e-13;
  const T rel_tol = 4e-4;
  const size_t ecc_size = 10;
  const size_t anom_size = 100;
  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    const kepler::starters::rppb starter(eccentricity);
    for (size_t m = 0; m < anom_size; ++m) {
      const T ecc_anom_expect = kepler::constants::pi<T>() * m / T(anom_size - 1);
      auto mean_anomaly = ecc_anom_expect - eccentricity * std::sin(ecc_anom_expect);
      auto ecc_anom = starter.start(mean_anomaly);
      if (m == 0) {
        REQUIRE_THAT(ecc_anom, WithinAbs(ecc_anom_expect, abs_tol));
      } else {
        REQUIRE_THAT(ecc_anom, WithinRel(ecc_anom_expect, rel_tol));
      }
    }
  }
}
