#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <string>

#include "kepler/kepler.hpp"
#include "reference/batman.hpp"
#include "reference/contour.hpp"
#include "reference/radvel.hpp"

template <typename R, typename S = kepler::starters::basic<typename R::value_type>>
struct Benchmark {
  typedef typename R::value_type value_type;
  typedef R refiner_type;
  typedef S starter_type;
};

template <typename S>
struct RefBenchmark {
  typedef typename S::value_type value_type;
  typedef S solver_type;
};

#define DEFAULT_NUM_DATA 1000
#define GENERATE_TEST_DATA(SIZE)                          \
  using T = typename TestType::value_type;                \
  std::vector<T> mean_anomaly(SIZE), ecc_anomaly(SIZE);   \
  for (size_t m = 0; m < SIZE; ++m) {                     \
    mean_anomaly[m] = T(100.) * m / T(SIZE - 1) - T(50.); \
  }

TEMPLATE_PRODUCT_TEST_CASE("Baseline", "[bench][baseline]", Benchmark,
                           (kepler::refiners::noop<float>, kepler::refiners::noop<double>)) {
  const size_t num_anom = DEFAULT_NUM_DATA;
  GENERATE_TEST_DATA(num_anom);
  BENCHMARK("Baseline") {
    for (std::size_t n = 0; n < num_anom; ++n) {
      ecc_anomaly[n] = std::sin(mean_anomaly[n]) + std::cos(mean_anomaly[n]);
    }
  };
}

#define MAIN_BENCHMARK(NAME, TAGS, ALGO)                                                        \
  TEMPLATE_PRODUCT_TEST_CASE(NAME, TAGS, Benchmark, (ALGO)) {                                   \
    const size_t num_ecc = 5;                                                                   \
    const size_t num_anom = DEFAULT_NUM_DATA;                                                   \
    const typename TestType::refiner_type refiner;                                              \
    for (size_t n = 0; n < num_ecc; ++n) {                                                      \
      GENERATE_TEST_DATA(num_anom);                                                             \
      const T eccentricity = n / T(num_ecc);                                                    \
      auto name = "e = " + std::to_string(eccentricity) + "; n = " + std::to_string(num_anom);  \
      BENCHMARK(name.c_str()) {                                                                 \
        return kepler::solve<typename TestType::starter_type, typename TestType::refiner_type>( \
            eccentricity, num_anom, mean_anomaly.data(), ecc_anomaly.data(), refiner);          \
      };                                                                                        \
    }                                                                                           \
  }

MAIN_BENCHMARK("First-order iterative float", "[bench][iterative][first-order][float]",
               (kepler::refiners::iterative<1, float>));
MAIN_BENCHMARK("First-order iterative double", "[bench][iterative][first-order][double]",
               (kepler::refiners::iterative<1, double>));

MAIN_BENCHMARK("Third-order iterative float", "[bench][iterative][third-order][float]",
               (kepler::refiners::iterative<3, float>));
MAIN_BENCHMARK("Third-order iterative double", "[bench][iterative][third-order][double]",
               (kepler::refiners::iterative<3, double>));

#undef MAIN_BENCHMARK

#define SIMD_BENCHMARK(NAME, TAGS, ALGO)                                                       \
  TEMPLATE_PRODUCT_TEST_CASE(NAME, TAGS, Benchmark, (ALGO)) {                                  \
    const size_t num_ecc = 5;                                                                  \
    const size_t num_anom = DEFAULT_NUM_DATA;                                                  \
    const typename TestType::refiner_type refiner;                                             \
    for (size_t n = 0; n < num_ecc; ++n) {                                                     \
      GENERATE_TEST_DATA(num_anom);                                                            \
      const T eccentricity = n / T(num_ecc);                                                   \
      auto name = "e = " + std::to_string(eccentricity) + "; n = " + std::to_string(num_anom); \
      BENCHMARK(name.c_str()) {                                                                \
        return kepler::solve_simd<typename TestType::starter_type,                             \
                                  typename TestType::refiner_type>(                            \
            eccentricity, num_anom, mean_anomaly.data(), ecc_anomaly.data(), refiner);         \
      };                                                                                       \
    }                                                                                          \
  }

SIMD_BENCHMARK("First-order iterative float (SIMD)",
               "[bench][iterative][first-order][float][simd]",
               (kepler::refiners::iterative<1, float>));
SIMD_BENCHMARK("First-order iterative double (SIMD)",
               "[bench][iterative][first-order][double][simd]",
               (kepler::refiners::iterative<1, double>));

SIMD_BENCHMARK("Third-order iterative float (SIMD)",
               "[bench][iterative][third-order][float][simd]",
               (kepler::refiners::iterative<3, float>));
SIMD_BENCHMARK("Third-order iterative double (SIMD)",
               "[bench][iterative][third-order][double][simd]",
               (kepler::refiners::iterative<3, double>));

#undef SIMD_BENCHMARK

#define REFERENCE_BENCHMARK(NAME, TAGS, ALGO)                                                  \
  TEMPLATE_PRODUCT_TEST_CASE(NAME, TAGS, RefBenchmark, (ALGO)) {                               \
    const size_t num_ecc = 5;                                                                  \
    const size_t num_anom = DEFAULT_NUM_DATA;                                                  \
    typename TestType::solver_type solver;                                                     \
    for (size_t n = 0; n < num_ecc; ++n) {                                                     \
      GENERATE_TEST_DATA(num_anom);                                                            \
      const T eccentricity = n / T(num_ecc);                                                   \
      auto name = "e = " + std::to_string(eccentricity) + "; n = " + std::to_string(num_anom); \
      BENCHMARK(name.c_str()) {                                                                \
        solver.setup(eccentricity);                                                            \
        for (std::size_t n = 0; n < num_anom; ++n) {                                           \
          ecc_anomaly[n] = solver.solve(mean_anomaly[n]);                                      \
        }                                                                                      \
        return ecc_anomaly[num_anom - 1];                                                      \
      };                                                                                       \
    }                                                                                          \
  }

REFERENCE_BENCHMARK("Batman", "[bench][reference][batman][first-order]",
                    kepler::reference::batman);
REFERENCE_BENCHMARK("RadVel", "[bench][reference][radvel][third-order]",
                    kepler::reference::radvel);
REFERENCE_BENCHMARK("Contour 8", "[bench][reference][contour][contour-8]",
                    kepler::reference::contour<8>);
REFERENCE_BENCHMARK("Contour 16", "[bench][reference][contour][contour-16]",
                    kepler::reference::contour<16>);

#undef REFERENCE_BENCHMARK
