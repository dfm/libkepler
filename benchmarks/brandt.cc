#include "kepler_benchmarks/brandt.h"

#include <iostream>

#include "kepler_benchmarks/benchmark.h"

int main() {
  kepler_benchmarks::Brandt<double> solver;
  kepler_benchmarks::run_general_benchmark("brandt.csv", solver);
  return 0;
}