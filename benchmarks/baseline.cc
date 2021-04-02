#include <iostream>

#define KEPLER_BENCHMARK_ECCENTRICITY_STEPS 1
#define KEPLER_BENCHMARK_ITERATIONS 10000000

#include "kepler_benchmarks/benchmark.h"
#include "kepler_benchmarks/iterative.h"

int main() {
  kepler_benchmarks::DummyBaselineSolver<double> baseline;
  kepler_benchmarks::run_general_benchmark("baseline.csv", baseline);
  return 0;
}