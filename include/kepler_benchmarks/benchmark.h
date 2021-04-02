#ifndef _KEPLER_BENCHMARKS_BENCHMARK_H_
#define _KEPLER_BENCHMARKS_BENCHMARK_H_

#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <vector>

#ifndef KEPLER_BENCHMARK_ECCENTRICITY_STEPS
#define KEPLER_BENCHMARK_ECCENTRICITY_STEPS 7
#endif

#ifndef KEPLER_BENCHMARK_ECCENTRIC_ANOMALY_STEPS
#define KEPLER_BENCHMARK_ECCENTRIC_ANOMALY_STEPS 13
#endif

#ifndef KEPLER_BENCHMARK_ITERATIONS
#define KEPLER_BENCHMARK_ITERATIONS 5000000
#endif

namespace kepler_benchmarks {

template <typename T>
struct BenchmarkResults {
  ssize_t x;
  ssize_t y;
  T eccentricity;
  T mean_anomaly;
  T true_eccentric_anomaly;
  T computed_sin_eccentric_anomaly;
  T computed_cos_eccentric_anomaly;
  ssize_t iterations;
  std::chrono::milliseconds elapsed;

  std::string to_csv() {
    std::ostringstream stream;
    stream << std::scientific;
    stream << x << ",";
    stream << y << ",";
    stream << eccentricity << ",";
    stream << eccentricity << ",";
    stream << eccentricity << ",";
    stream << mean_anomaly << ",";
    stream << true_eccentric_anomaly << ",";
    stream << computed_sin_eccentric_anomaly << ",";
    stream << computed_cos_eccentric_anomaly << ",";
    stream << fabs(sin(true_eccentric_anomaly) - computed_sin_eccentric_anomaly) << ",";
    stream << fabs(cos(true_eccentric_anomaly) - computed_cos_eccentric_anomaly) << ",";
    stream << iterations << ",";
    stream << elapsed.count() << "\n";
    return stream.str();
  }

  static std::string header() {
    std::ostringstream stream;
    stream << "index_ecc,index_anom,";
    stream << "eccentricity,mean_anomaly,true_eccentric_anomaly,";
    stream << "computed_sin,computed_cos,error_sin,error_cos,";
    stream << "iterations,time_ms\n";
    return stream.str();
  }
};

template <typename SolverType>
void run_general_benchmark(const std::string& filename, SolverType& solver) {
  typedef typename SolverType::Scalar T;

  const ssize_t eccentricity_steps = KEPLER_BENCHMARK_ECCENTRICITY_STEPS;
  const ssize_t eccentric_anomaly_steps = KEPLER_BENCHMARK_ECCENTRIC_ANOMALY_STEPS;
  const ssize_t iterations = KEPLER_BENCHMARK_ITERATIONS;

  std::vector<BenchmarkResults<T>> results;
  T sinE = 0, cosE = 0;

  std::ofstream output_file_stream;
  output_file_stream.open(filename);
  output_file_stream << BenchmarkResults<T>::header();

  for (ssize_t i = 0; i < eccentricity_steps; ++i) {
    T eccentricity = T(i) / T(eccentricity_steps);
    for (ssize_t j = 0; j < eccentric_anomaly_steps; ++j) {
      T eccentric_anomaly = 2 * M_PI * T(j) / T(eccentric_anomaly_steps);

      // Compute the mean anomaly from the true eccentric anomaly and wrap into [0, 2*pi)
      T mean_anomaly = eccentric_anomaly - eccentricity * sin(eccentric_anomaly);
      while (mean_anomaly >= 2 * M_PI) mean_anomaly -= 2 * M_PI;
      while (mean_anomaly < 0) mean_anomaly += 2 * M_PI;

      // Time several loops
      auto start = std::chrono::high_resolution_clock::now();
      for (ssize_t k = 0; k < iterations; ++k) {
        solver.evaluate(mean_anomaly, eccentricity, &sinE, &cosE);
      }
      auto end = std::chrono::high_resolution_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

      // Save the results
      BenchmarkResults<T> results(
          {i, j, eccentricity, mean_anomaly, eccentric_anomaly, sinE, cosE, iterations, elapsed});
      output_file_stream << results.to_csv();
    }
  }
  output_file_stream.close();
}

template <typename SolverType>
void run_fixed_eccentricity_benchmark(const std::string& filename, SolverType& solver) {
  typedef typename SolverType::Scalar T;

  const ssize_t eccentricity_steps = KEPLER_BENCHMARK_ECCENTRICITY_STEPS;
  const ssize_t eccentric_anomaly_steps = KEPLER_BENCHMARK_ECCENTRIC_ANOMALY_STEPS;
  const ssize_t iterations = KEPLER_BENCHMARK_ITERATIONS;

  std::vector<BenchmarkResults<T>> results;
  T sinE = 0, cosE = 0;

  std::ofstream output_file_stream;
  output_file_stream.open(filename);
  output_file_stream << BenchmarkResults<T>::header();

  for (ssize_t i = 0; i < eccentricity_steps; ++i) {
    T eccentricity = T(i) / T(eccentricity_steps);
    solver.precompute_for_eccentricity(eccentricity);

    for (ssize_t j = 0; j < eccentric_anomaly_steps; ++j) {
      T eccentric_anomaly = 2 * M_PI * T(j) / T(eccentric_anomaly_steps);

      // Compute the mean anomaly from the true eccentric anomaly and wrap into [0, 2*pi)
      T mean_anomaly = eccentric_anomaly - eccentricity * sin(eccentric_anomaly);
      while (mean_anomaly >= 2 * M_PI) mean_anomaly -= 2 * M_PI;
      while (mean_anomaly < 0) mean_anomaly += 2 * M_PI;

      // Time several loops
      auto start = std::chrono::high_resolution_clock::now();
      for (ssize_t k = 0; k < iterations; ++k) {
        solver.evaluate(mean_anomaly, &sinE, &cosE);
      }
      auto end = std::chrono::high_resolution_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

      // Save the results
      BenchmarkResults<T> results(
          {i, j, eccentricity, mean_anomaly, eccentric_anomaly, sinE, cosE, iterations, elapsed});
      output_file_stream << results.to_csv();
    }
  }
  output_file_stream.close();
}

}  // namespace kepler_benchmarks

#endif
