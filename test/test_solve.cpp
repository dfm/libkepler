#include <cmath>
#include <vector>

#include "./test_utils.hpp"
#include "kepler/kepler.hpp"

TEMPLATE_PRODUCT_TEST_CASE(
    "SIMD comparison", "[refiners][simd]", SolveTestCase,
    (kepler::refiners::noop<double>, (kepler::refiners::iterative<1, float>),
     (kepler::refiners::iterative<1, double>), (kepler::refiners::iterative<2, double>),
     (kepler::refiners::iterative<3, double>), (kepler::refiners::iterative<4, double>),
     (kepler::refiners::iterative<5, double>), (kepler::refiners::iterative<6, double>),
     (kepler::refiners::iterative<7, double>),
     (kepler::refiners::non_iterative<4, double>, kepler::starters::markley<double>),
     (kepler::refiners::non_iterative<4, float>, kepler::starters::markley<float>))) {
  using T = typename TestType::value_type;
  const T abs_tol = tolerance<TestType>::abs;
  const size_t ecc_size = 10;
  const size_t anom_size = 1000;
  const typename TestType::refiner_type refiner;
  std::vector<T> mean_anomaly(anom_size), ecc_anom(anom_size), ecc_anom_simd(anom_size);

  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    for (size_t m = 0; m < anom_size; ++m) {
      mean_anomaly[m] = 100. * m / T(anom_size - 1) - 50.;
    }

    kepler::solve<typename TestType::starter_type, typename TestType::refiner_type>(
        eccentricity, anom_size, mean_anomaly.data(), ecc_anom.data(), refiner);
    kepler::solve_simd<typename TestType::starter_type, typename TestType::refiner_type>(
        eccentricity, anom_size, mean_anomaly.data(), ecc_anom_simd.data(), refiner);

    for (size_t m = 0; m < anom_size; ++m) {
      REQUIRE_THAT(ecc_anom_simd[m], WithinAbs(ecc_anom[m], abs_tol));
    }
  }
}
