#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <string>

#include "kepler/kepler.hpp"

template <typename R, typename S = kepler::starters::basic<typename R::value_type>>
struct Benchmark {
  typedef typename R::value_type value_type;
  typedef R refiner_type;
  typedef S starter_type;
};

#define DEFAULT_NUM_DATA 1000
#define GENERATE_TEST_DATA(SIZE)                          \
  using T = typename TestType::value_type;                \
  std::vector<T> mean_anomaly(SIZE), ecc_anomaly(SIZE);   \
  for (size_t m = 0; m < SIZE; ++m) {                     \
    mean_anomaly[m] = T(100.) * m / T(SIZE - 1) - T(50.); \
  }

TEMPLATE_PRODUCT_TEST_CASE("Baseline", "[baseline][bench]", Benchmark,
                           (kepler::refiners::noop<float>, kepler::refiners::noop<double>)) {
  const size_t num_anom = DEFAULT_NUM_DATA;
  GENERATE_TEST_DATA(num_anom);
  BENCHMARK("Baseline") {
    for (std::size_t n = 0; n < num_anom; ++n) {
      ecc_anomaly[n] = std::sin(mean_anomaly[n]) + std::cos(mean_anomaly[n]);
    }
  };
}

TEMPLATE_PRODUCT_TEST_CASE("Iterative", "[iterative][bench]", Benchmark,
                           ((kepler::refiners::iterative<1, float>),
                            (kepler::refiners::iterative<1, double>))) {
  const size_t num_ecc = 5;
  const size_t num_anom = DEFAULT_NUM_DATA;
  const typename TestType::refiner_type refiner;
  for (size_t n = 0; n < num_ecc; ++n) {
    GENERATE_TEST_DATA(num_anom);
    const T eccentricity = n / T(num_ecc);
    auto name = "Eccentricity = " + std::to_string(eccentricity);
    BENCHMARK(name.c_str()) {
      return kepler::solve<typename TestType::starter_type, typename TestType::refiner_type>(
          eccentricity, num_anom, mean_anomaly.data(), ecc_anomaly.data(), refiner);
    };
  }
}

TEMPLATE_PRODUCT_TEST_CASE("Iterative (SIMD)", "[iterative][bench][simd]", Benchmark,
                           ((kepler::refiners::iterative<1, float>),
                            (kepler::refiners::iterative<1, double>))) {
  const size_t num_ecc = 5;
  const size_t num_anom = DEFAULT_NUM_DATA;
  const typename TestType::refiner_type refiner;
  for (size_t n = 0; n < num_ecc; ++n) {
    GENERATE_TEST_DATA(num_anom);
    const T eccentricity = n / T(num_ecc);
    auto name = "Eccentricity = " + std::to_string(eccentricity);
    BENCHMARK(name.c_str()) {
      return kepler::solve_simd<typename TestType::starter_type, typename TestType::refiner_type>(
          eccentricity, num_anom, mean_anomaly.data(), ecc_anomaly.data(), refiner);
    };
  }
}
