#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iomanip>
#include <sstream>

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
#define GENERATE_TEST_DATA(SIZE)                                                                \
  using T = typename TestType::value_type;                                                      \
  std::vector<T> mean_anomaly(SIZE), ecc_anomaly(SIZE), sin_ecc_anom(SIZE), cos_ecc_anom(SIZE); \
  for (size_t m = 0; m < SIZE; ++m) {                                                           \
    mean_anomaly[m] = T(100.) * m / T(SIZE - 1) - T(50.);                                       \
  }

TEMPLATE_PRODUCT_TEST_CASE("baselinef", "[bench][baseline][float]", Benchmark,
                           kepler::refiners::noop<float>) {
  const size_t num_anom = DEFAULT_NUM_DATA;
  GENERATE_TEST_DATA(num_anom);
  BENCHMARK("e=0; n=1000") {
    for (std::size_t n = 0; n < num_anom; ++n) {
      ecc_anomaly[n] = std::sin(mean_anomaly[n]) + std::cos(mean_anomaly[n]);
    }
  };
}

TEMPLATE_PRODUCT_TEST_CASE("baselined", "[bench][baseline][double]", Benchmark,
                           kepler::refiners::noop<double>) {
  const size_t num_anom = DEFAULT_NUM_DATA;
  GENERATE_TEST_DATA(num_anom);
  BENCHMARK("e=0; n=1000") {
    for (std::size_t n = 0; n < num_anom; ++n) {
      ecc_anomaly[n] = std::sin(mean_anomaly[n]) + std::cos(mean_anomaly[n]);
    }
  };
}

#define MAIN_BENCHMARK(NAME, TAGS, ALGO)                                                          \
  TEMPLATE_PRODUCT_TEST_CASE(NAME, TAGS, Benchmark, (ALGO)) {                                     \
    const size_t num_ecc = 5;                                                                     \
    const size_t num_anom = DEFAULT_NUM_DATA;                                                     \
    const typename TestType::refiner_type refiner;                                                \
    for (size_t n = 0; n < num_ecc; ++n) {                                                        \
      GENERATE_TEST_DATA(num_anom);                                                               \
      const T eccentricity = (T(n) + T(0.5)) / T(num_ecc);                                        \
      std::ostringstream name;                                                                    \
      name << std::setprecision(1) << "e=" << eccentricity << "; n=" << num_anom;                 \
      BENCHMARK(name.str().c_str()) {                                                             \
        return kepler::solver::solve<typename TestType::starter_type,                             \
                                     typename TestType::refiner_type>(                            \
            eccentricity, num_anom, mean_anomaly.data(), ecc_anomaly.data(), sin_ecc_anom.data(), \
            cos_ecc_anom.data(), refiner);                                                        \
      };                                                                                          \
    }                                                                                             \
  }

MAIN_BENCHMARK("iter1f", "[bench][iterative][first-order][float]",
               (kepler::refiners::iterative<1, float>))
MAIN_BENCHMARK("iter1d", "[bench][iterative][first-order][double]",
               (kepler::refiners::iterative<1, double>))

MAIN_BENCHMARK("iter3f", "[bench][iterative][third-order][float]",
               (kepler::refiners::iterative<3, float>))
MAIN_BENCHMARK("iter3d", "[bench][iterative][third-order][double]",
               (kepler::refiners::iterative<3, double>))

MAIN_BENCHMARK("markley95f", "[bench][non-iterative][markley][float]",
               (kepler::refiners::non_iterative<3, float>, kepler::starters::markley<float>))
MAIN_BENCHMARK("markley95d", "[bench][non-iterative][markley][double]",
               (kepler::refiners::non_iterative<3, double>, kepler::starters::markley<double>))

MAIN_BENCHMARK("brandt21f", "[bench][non-iterative][brandt][float]",
               (kepler::refiners::brandt<float>, kepler::starters::raposo_pulido_brandt<float>))
MAIN_BENCHMARK("brandt21d", "[bench][non-iterative][brandt][double]",
               (kepler::refiners::brandt<double>, kepler::starters::raposo_pulido_brandt<double>))

#undef MAIN_BENCHMARK

#define SIMD_BENCHMARK(NAME, TAGS, ALGO)                                                          \
  TEMPLATE_PRODUCT_TEST_CASE(NAME, TAGS, Benchmark, (ALGO)) {                                     \
    const size_t num_ecc = 5;                                                                     \
    const size_t num_anom = DEFAULT_NUM_DATA;                                                     \
    const typename TestType::refiner_type refiner;                                                \
    for (size_t n = 0; n < num_ecc; ++n) {                                                        \
      GENERATE_TEST_DATA(num_anom);                                                               \
      const T eccentricity = (T(n) + T(0.5)) / T(num_ecc);                                        \
      std::ostringstream name;                                                                    \
      name << std::setprecision(1) << "e=" << eccentricity << "; n=" << num_anom;                 \
      BENCHMARK(name.str().c_str()) {                                                             \
        return kepler::solver::solve_simd<typename TestType::starter_type,                        \
                                          typename TestType::refiner_type>(                       \
            eccentricity, num_anom, mean_anomaly.data(), ecc_anomaly.data(), sin_ecc_anom.data(), \
            cos_ecc_anom.data(), refiner);                                                        \
      };                                                                                          \
    }                                                                                             \
  }

SIMD_BENCHMARK("iter1fv", "[bench][iterative][first-order][float][simd]",
               (kepler::refiners::iterative<1, float>))
SIMD_BENCHMARK("iter1dv", "[bench][iterative][first-order][double][simd]",
               (kepler::refiners::iterative<1, double>))

SIMD_BENCHMARK("iter3fv", "[bench][iterative][third-order][float][simd]",
               (kepler::refiners::iterative<3, float>))
SIMD_BENCHMARK("iter3dv", "[bench][iterative][third-order][double][simd]",
               (kepler::refiners::iterative<3, double>))

SIMD_BENCHMARK("markley95fv", "[bench][non-iterative][markley][float]",
               (kepler::refiners::non_iterative<3, float>, kepler::starters::markley<float>))
SIMD_BENCHMARK("markley95dv", "[bench][non-iterative][markley][double]",
               (kepler::refiners::non_iterative<3, double>, kepler::starters::markley<double>))

SIMD_BENCHMARK("brandt21fv", "[bench][non-iterative][brandt][float][simd]",
               (kepler::refiners::brandt<float>, kepler::starters::raposo_pulido_brandt<float>))
SIMD_BENCHMARK("brandt21dv", "[bench][non-iterative][brandt][double][simd]",
               (kepler::refiners::brandt<double>, kepler::starters::raposo_pulido_brandt<double>))

#undef SIMD_BENCHMARK

#define REFERENCE_BENCHMARK(NAME, TAGS, ALGO)                                         \
  TEMPLATE_PRODUCT_TEST_CASE(NAME, TAGS, RefBenchmark, (ALGO)) {                      \
    const size_t num_ecc = 5;                                                         \
    const size_t num_anom = DEFAULT_NUM_DATA;                                         \
    typename TestType::solver_type solver;                                            \
    for (size_t n = 0; n < num_ecc; ++n) {                                            \
      GENERATE_TEST_DATA(num_anom);                                                   \
      const T eccentricity = (T(n) + T(0.5)) / T(num_ecc);                            \
      std::ostringstream name;                                                        \
      name << std::setprecision(1) << "e=" << eccentricity << "; n=" << num_anom;     \
      BENCHMARK(name.str().c_str()) {                                                 \
        solver.setup(eccentricity);                                                   \
        for (std::size_t n = 0; n < num_anom; ++n) {                                  \
          ecc_anomaly[n] = solver.solve(mean_anomaly[n]);                             \
          sin_ecc_anom[n] = std::sin(ecc_anomaly[n]);                                 \
          cos_ecc_anom[n] = std::cos(ecc_anomaly[n]);                                 \
        }                                                                             \
        return std::make_tuple(ecc_anomaly[num_anom - 1], sin_ecc_anom[num_anom - 1], \
                               cos_ecc_anom[num_anom - 1]);                           \
      };                                                                              \
    }                                                                                 \
  }

REFERENCE_BENCHMARK("ref:batman", "[bench][reference][batman][iterative][first-order][double]",
                    kepler::reference::batman)
REFERENCE_BENCHMARK("ref:radvel", "[bench][reference][radvel][iterative][third-order][double]",
                    kepler::reference::radvel)
REFERENCE_BENCHMARK("ref:contour8",
                    "[bench][reference][non-iterative][contour][contour-8][double]",
                    kepler::reference::contour<8>)
REFERENCE_BENCHMARK("ref:contour16",
                    "[bench][reference][non-iterative][contour][contour-16][double]",
                    kepler::reference::contour<16>)

#undef REFERENCE_BENCHMARK
