#include <cmath>
#include <vector>

#include "./test_utils.hpp"
#include "kepler/kepler/constants.hpp"
#include "kepler/kepler/refiners.hpp"
#include "kepler/kepler/solver.hpp"
#include "kepler/kepler/starters.hpp"

using namespace kepler;
template <typename T>
struct has_tolerance {
  static const bool value = false;
};

template <int order, typename T>
struct has_tolerance<refiners::iterative<order, T>> {
  static const bool value = true;
};

TEMPLATE_PRODUCT_TEST_CASE("Refiners", "[refiners]", SolveTestCase,
                           ((refiners::iterative<1, float>), (refiners::iterative<1, double>),
                            (refiners::iterative<2, double>), (refiners::iterative<3, double>),
                            (refiners::iterative<4, double>), (refiners::iterative<5, double>),
                            (refiners::iterative<6, double>), (refiners::iterative<7, double>),
                            (refiners::non_iterative<3, double>, starters::markley<double>),
                            (refiners::non_iterative<3, float>, starters::markley<float>),
                            (refiners::brandt<float>, starters::raposo_pulido_brandt<float>),
                            (refiners::brandt<double>, starters::raposo_pulido_brandt<double>))) {
  using T = typename TestType::value_type;
  const size_t ecc_size = 10;
  const size_t anom_size = 1000;
  const typename TestType::refiner_type refiner;
  std::vector<T> ecc_anom_expect(anom_size), mean_anomaly(anom_size), ecc_anom_calc(anom_size),
      sin_ecc_anom(anom_size), cos_ecc_anom(anom_size);

  T abs_tol = tolerance<TestType>::abs;
  if constexpr (has_tolerance<typename TestType::refiner_type>::value) {
    abs_tol = 20 * refiner.tolerance;
  }

  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    for (size_t m = 0; m < anom_size; ++m) {
      ecc_anom_expect[m] = T(100.) * m / T(anom_size - 1) - T(50.);
      mean_anomaly[m] = ecc_anom_expect[m] - eccentricity * std::sin(ecc_anom_expect[m]);
    }

    solver::solve<typename TestType::starter_type, typename TestType::refiner_type>(
        eccentricity, anom_size, mean_anomaly.data(), ecc_anom_calc.data(), sin_ecc_anom.data(),
        cos_ecc_anom.data(), refiner);

    for (size_t m = 0; m < anom_size; ++m) {
      REQUIRE_THAT(std::sin(ecc_anom_calc[m]), WithinAbs(std::sin(ecc_anom_expect[m]), abs_tol));
      REQUIRE_THAT(std::cos(ecc_anom_calc[m]), WithinAbs(std::cos(ecc_anom_expect[m]), abs_tol));
      REQUIRE_THAT(sin_ecc_anom[m], WithinAbs(std::sin(ecc_anom_expect[m]), abs_tol));
      REQUIRE_THAT(cos_ecc_anom[m], WithinAbs(std::cos(ecc_anom_expect[m]), abs_tol));
    }
  }
}

TEMPLATE_PRODUCT_TEST_CASE("SIMD comparison", "[refiners][simd]", SolveTestCase,
                           (refiners::noop<double>, (refiners::iterative<1, float>),
                            (refiners::iterative<1, double>), (refiners::iterative<2, double>),
                            (refiners::iterative<3, double>), (refiners::iterative<4, double>),
                            (refiners::iterative<5, double>), (refiners::iterative<6, double>),
                            (refiners::iterative<7, double>),
                            (refiners::non_iterative<3, double>, starters::markley<double>),
                            (refiners::non_iterative<3, float>, starters::markley<float>),
                            (refiners::brandt<float>, starters::raposo_pulido_brandt<float>),
                            (refiners::brandt<double>, starters::raposo_pulido_brandt<double>))) {
  using T = typename TestType::value_type;
  using B = xs::batch<T>;
  constexpr std::size_t simd_size = B::size;
  const T abs_tol = tolerance<TestType>::abs;
  const size_t ecc_size = 10;
  const size_t anom_size = 100 * simd_size;
  alignas(B::arch_type::alignment()) std::array<T, simd_size> ecc_anom, mean_anom;

  const typename TestType::refiner_type refiner;
  std::vector<T> ecc_anom_expect(anom_size), mean_anomaly(anom_size), ecc_anom_calc(anom_size);

  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    const typename TestType::starter_type starter(eccentricity);

    for (size_t m = 0; m < anom_size; m += simd_size) {
      for (size_t k = 0; k < simd_size; ++k) {
        mean_anom[k] = constants::pi<T>() * (m + k) / T(anom_size - 1);
      }
      auto mean_anom_b = xs::load_aligned(mean_anom.data());
      auto ecc_anom_b = starter.start(mean_anom_b);
      ecc_anom_b = refiner.refine(eccentricity, mean_anom_b, ecc_anom_b);
      ecc_anom_b.store_aligned(ecc_anom.data());

      for (size_t k = 0; k < simd_size; ++k) {
        auto expect = refiner.refine(eccentricity, mean_anom[k], starter.start(mean_anom[k]));
        REQUIRE_THAT(ecc_anom[k], WithinAbs(expect, abs_tol));
      }
    }
  }
}
