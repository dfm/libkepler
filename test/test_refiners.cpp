#include <cmath>
#include <vector>

#include "./test_utils.hpp"
#include "kepler/kepler.hpp"

TEMPLATE_PRODUCT_TEST_CASE(
    "Refiners", "[refiners]", SolveTestCase,
    ((kepler::refiners::iterative<1, float>), (kepler::refiners::iterative<1, double>),
     (kepler::refiners::iterative<2, double>), (kepler::refiners::iterative<3, double>),
     (kepler::refiners::iterative<4, double>), (kepler::refiners::iterative<5, double>),
     (kepler::refiners::iterative<6, double>), (kepler::refiners::iterative<7, double>),
     (kepler::refiners::non_iterative<4, double>, kepler::starters::markley<double>),
     (kepler::refiners::non_iterative<4, float>, kepler::starters::markley<float>))) {
  using T = typename TestType::value_type;
  const T abs_tol = tolerance<TestType>::abs;
  const size_t ecc_size = 10;
  const size_t anom_size = 1000;
  const typename TestType::refiner_type refiner;
  std::vector<T> ecc_anom_expect(anom_size), mean_anomaly(anom_size), ecc_anom_calc(anom_size);

  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    for (size_t m = 0; m < anom_size; ++m) {
      ecc_anom_expect[m] = 100. * m / T(anom_size - 1) - 50.;
      mean_anomaly[m] = ecc_anom_expect[m] - eccentricity * std::sin(ecc_anom_expect[m]);
    }

    kepler::solve<typename TestType::starter_type, typename TestType::refiner_type>(
        eccentricity, anom_size, mean_anomaly.data(), ecc_anom_calc.data(), refiner);

    for (size_t m = 0; m < anom_size; ++m) {
      REQUIRE_THAT(std::sin(ecc_anom_calc[m]), WithinAbs(std::sin(ecc_anom_expect[m]), abs_tol));
      REQUIRE_THAT(std::cos(ecc_anom_calc[m]), WithinAbs(std::cos(ecc_anom_expect[m]), abs_tol));
    }
  }
}

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
  constexpr std::size_t simd_size = xs::simd_type<T>::size;
  const T abs_tol = tolerance<TestType>::abs;
  const size_t ecc_size = 10;
  const size_t anom_size = 100 * simd_size;
  alignas(xs::batch<T>::arch_type::alignment()) std::array<T, simd_size> ecc_anom, mean_anom;

  const typename TestType::refiner_type refiner;
  std::vector<T> ecc_anom_expect(anom_size), mean_anomaly(anom_size), ecc_anom_calc(anom_size);

  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    const typename TestType::starter_type starter(eccentricity);

    for (size_t m = 0; m < anom_size; m += simd_size) {
      for (size_t k = 0; k < simd_size; ++k) {
        mean_anom[k] = kepler::constants::pi<T>() * (m + k) / T(anom_size - 1);
      }
      auto mean_anom_b = xs::load_aligned(mean_anom.data());
      auto ecc_anom_b = starter.start(mean_anom_b);
      ecc_anom_b = refiner.refine(eccentricity, mean_anom_b, ecc_anom_b);
      ecc_anom_b.store_aligned(ecc_anom.data());

      for (size_t k = 0; k < simd_size; ++k) {
        REQUIRE_THAT(ecc_anom[k], WithinAbs(refiner.refine(eccentricity, mean_anom[k],
                                                           starter.start(mean_anom[k])),
                                            abs_tol));
      }
    }
  }
}
