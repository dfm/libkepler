#include <iostream>

#include "kepler_benchmarks/benchmark.h"
#include "kepler_benchmarks/brandt.h"

int main() {
  kepler_benchmarks::BrandtFixed<double> solver;
  kepler_benchmarks::run_fixed_eccentricity_benchmark("brandt_fixed.csv", solver);
  return 0;
}