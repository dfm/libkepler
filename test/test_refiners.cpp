#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <vector>

#include "kepler/kepler.hpp"

using namespace Catch::Matchers;

TEMPLATE_TEST_CASE("Iterative refiners", "[refiners]", (kepler::refiners::iterative<1, double>),
                   (kepler::refiners::iterative<2, double>),
                   (kepler::refiners::iterative<3, double>),
                   (kepler::refiners::iterative<4, double>),
                   (kepler::refiners::iterative<5, double>),
                   (kepler::refiners::iterative<6, double>),
                   (kepler::refiners::iterative<7, double>)) {
  using T = double;
  const T abs_tol = 5e-12;
  const size_t ecc_size = 10;
  const size_t anom_size = 1000;
  const TestType refiner;
  std::vector<T> ecc_anom_expect(anom_size), mean_anomaly(anom_size), ecc_anom_calc(anom_size);

  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    for (size_t m = 0; m < anom_size; ++m) {
      ecc_anom_expect[m] = 100. * m / T(anom_size - 1) - 50.;
      mean_anomaly[m] = ecc_anom_expect[m] - eccentricity * std::sin(ecc_anom_expect[m]);
    }

    kepler::solve(eccentricity, anom_size, mean_anomaly.data(), ecc_anom_calc.data(), refiner);

    for (size_t m = 0; m < anom_size; ++m) {
      REQUIRE_THAT(std::sin(ecc_anom_calc[m]), WithinAbs(std::sin(ecc_anom_expect[m]), abs_tol));
      REQUIRE_THAT(std::cos(ecc_anom_calc[m]), WithinAbs(std::cos(ecc_anom_expect[m]), abs_tol));
    }
  }
}

TEST_CASE("Nijenhuis / Markley", "[refiners]") {
  using T = double;
  const T abs_tol = 5e-12;
  const size_t ecc_size = 10;
  const size_t anom_size = 1000;
  std::vector<T> ecc_anom_expect(anom_size), mean_anomaly(anom_size), ecc_anom_calc(anom_size);

  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    for (size_t m = 0; m < anom_size; ++m) {
      ecc_anom_expect[m] = 100. * m / T(anom_size - 1) - 50.;
      mean_anomaly[m] = ecc_anom_expect[m] - eccentricity * std::sin(ecc_anom_expect[m]);
    }

    kepler::solve_non_iterative_1(eccentricity, anom_size, mean_anomaly.data(),
                                  ecc_anom_calc.data());

    for (size_t m = 0; m < anom_size; ++m) {
      REQUIRE_THAT(std::sin(ecc_anom_calc[m]), WithinAbs(std::sin(ecc_anom_expect[m]), abs_tol));
      REQUIRE_THAT(std::cos(ecc_anom_calc[m]), WithinAbs(std::cos(ecc_anom_expect[m]), abs_tol));
    }
  }
}
