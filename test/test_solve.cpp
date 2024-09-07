#include <cmath>
#include <vector>

#include "./test_utils.hpp"
#include "kepler/kepler/constants.hpp"
#include "kepler/kepler/refiners.hpp"
#include "kepler/kepler/solver.hpp"
#include "kepler/kepler/starters.hpp"

using namespace kepler;

TEMPLATE_PRODUCT_TEST_CASE("SIMD comparison", "[refiners][simd]", SolveTestCase,
                           ((refiners::noop<double>), (refiners::iterative<1, float>),
                            (refiners::iterative<1, double>), (refiners::iterative<2, double>),
                            (refiners::iterative<3, double>), (refiners::iterative<4, double>),
                            (refiners::iterative<5, double>), (refiners::iterative<6, double>),
                            (refiners::iterative<7, double>),
                            (refiners::non_iterative<3, double>, starters::markley<double>),
                            (refiners::non_iterative<3, float>, starters::markley<float>),
                            (refiners::brandt<float>, starters::raposo_pulido_brandt<float>),
                            (refiners::brandt<double>, starters::raposo_pulido_brandt<double>))) {
  using T = typename TestType::value_type;
  const T abs_tol = tolerance<TestType>::abs;
  const size_t ecc_size = 10;
  const size_t anom_size = 1003;
  const typename TestType::refiner_type refiner;
  std::vector<T> mean_anomaly(anom_size), ecc_anom(anom_size), sin_ecc_anom(anom_size),
      cos_ecc_anom(anom_size), ecc_anom_simd(anom_size), sin_ecc_anom_simd(anom_size),
      cos_ecc_anom_simd(anom_size);
  for (size_t m = 0; m < anom_size; ++m) {
    mean_anomaly[m] = T(100.) * m / T(anom_size - 1) - T(50.);
  }

  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);

    solver::solve<typename TestType::starter_type, typename TestType::refiner_type>(
        eccentricity, anom_size, mean_anomaly.data(), ecc_anom.data(), sin_ecc_anom.data(),
        cos_ecc_anom.data(), refiner);
    solver::solve_simd<typename TestType::starter_type, typename TestType::refiner_type>(
        eccentricity, anom_size, mean_anomaly.data(), ecc_anom_simd.data(),
        sin_ecc_anom_simd.data(), cos_ecc_anom_simd.data(), refiner);

    for (size_t m = 0; m < anom_size; ++m) {
      REQUIRE_THAT(ecc_anom_simd[m], WithinAbs(ecc_anom[m], abs_tol));
      REQUIRE_THAT(sin_ecc_anom_simd[m], WithinAbs(sin_ecc_anom[m], abs_tol));
      REQUIRE_THAT(cos_ecc_anom_simd[m], WithinAbs(cos_ecc_anom[m], abs_tol));
    }
  }
}
