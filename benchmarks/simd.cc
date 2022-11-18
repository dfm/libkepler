#include "kepler_benchmarks/simd.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include "kepler_benchmarks/benchmark.h"
#include "xsimd/xsimd.hpp"

namespace xs = xsimd;

int main() {
  using kepler_benchmarks::BenchmarkResults;
  using T = double;
  using vector_type = std::vector<T, xs::aligned_allocator<T>>;

  const std::size_t ecc_steps = 15;
  const std::size_t iterations = 500;
  const std::size_t size = 10000;

  std::vector<BenchmarkResults<T>> results;
  std::ofstream output_file_stream;
  output_file_stream.open("simd.csv");
  output_file_stream << BenchmarkResults<T>::header();

  vector_type ecc_anom_expect(size), ecc_anom_calc(size), mean_anom(size);
  for (std::size_t ei = 0; ei < ecc_steps; ++ei) {
    T ecc = T(ei + 0.5) / T(ecc_steps);

    for (std::size_t n = 0; n < size; ++n) {
      ecc_anom_expect[n] = T(2 * M_PI * (n + 0.5)) / T(size);
      mean_anom[n] = ecc_anom_expect[n] - ecc * std::sin(ecc_anom_expect[n]);
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t k = 0; k < iterations; ++k) {
      kepler::simd::solve<T, xs::aligned_mode>(size, ecc, mean_anom.data(), ecc_anom_calc.data());
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    auto max_diff = 0.0;
    auto max_diff_idx = 0;
    for (std::size_t n = 0; n < size; ++n) {
      auto diff = std::abs(std::sin(ecc_anom_calc[n]) - std::sin(ecc_anom_expect[n]));
      if (diff > max_diff) {
        max_diff = diff;
        max_diff_idx = n;
      }
    }
    std::cout << "ecc = " << ecc << ", max diff = " << max_diff << ", idx = " << max_diff_idx
              << std::endl;

    BenchmarkResults<T> results({static_cast<ssize_t>(ei), 0, ecc, T(0.0), T(0.0), T(0.0), T(0.0),
                                 iterations * size, elapsed});
    output_file_stream << results.to_csv();
  }
  output_file_stream.close();

  return 0;
}