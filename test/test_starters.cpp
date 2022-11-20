#include <array>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <xsimd/xsimd.hpp>

#include "kepler/kepler.hpp"

using namespace Catch::Matchers;
namespace xs = xsimd;

template <typename T>
struct default_abs {
  constexpr static T value = 1e-13;
};

template <>
struct default_abs<float> {
  constexpr static float value = 1e-5;
};

template <typename T>
struct default_rel {
  constexpr static T value = 5e-4;
};

template <>
struct default_rel<float> {
  constexpr static float value = 5e-2;
};

template <typename T>
struct tolerance {
  constexpr static typename T::value_type abs = default_abs<typename T::value_type>::value;
  constexpr static typename T::value_type rel = default_rel<typename T::value_type>::value;
};

template <>
struct tolerance<kepler::starters::mikkola<double>> {
  constexpr static double abs = 1e-13;
  constexpr static double rel = 0.002;
};

TEMPLATE_PRODUCT_TEST_CASE("Starters", "[starters]",
                           (kepler::starters::mikkola, kepler::starters::markley,
                            kepler::starters::rppb),
                           (double, float)) {
  using T = typename TestType::value_type;
  const T abs_tol = tolerance<TestType>::abs;
  const T rel_tol = tolerance<TestType>::rel;
  const size_t ecc_size = 10;
  const size_t anom_size = 100;
  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    const TestType starter(eccentricity);
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

TEMPLATE_PRODUCT_TEST_CASE("SIMD comparison", "[starters][simd]",
                           (kepler::starters::noop, kepler::starters::basic,
                            kepler::starters::mikkola, kepler::starters::markley),
                           (double, float)) {
  using T = typename TestType::value_type;
  constexpr std::size_t simd_size = xs::simd_type<T>::size;
  const T abs_tol = default_abs<T>::value;
  const size_t ecc_size = 10;
  const size_t anom_size = 100 * simd_size;
  alignas(xs::batch<T>::arch_type::alignment()) std::array<T, simd_size> ecc_anom, mean_anom;

  for (size_t n = 0; n < ecc_size; ++n) {
    const T eccentricity = n / T(ecc_size);
    const TestType starter(eccentricity);

    for (size_t m = 0; m < anom_size; m += simd_size) {
      for (size_t k = 0; k < simd_size; ++k) {
        mean_anom[k] = kepler::constants::pi<T>() * (m + k) / T(anom_size - 1);
      }
      auto mean_anom_b = xs::load_aligned(mean_anom.data());
      auto ecc_anom_b = starter.start(mean_anom_b);
      ecc_anom_b.store_aligned(ecc_anom.data());

      for (size_t k = 0; k < simd_size; ++k) {
        REQUIRE_THAT(ecc_anom[k], WithinAbs(starter.start(mean_anom[k]), abs_tol));
      }
    }
  }
}
