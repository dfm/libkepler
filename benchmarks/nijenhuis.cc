#include "kepler_benchmarks/nijenhuis.h"

#include <iostream>

#include "kepler_benchmarks/benchmark.h"

int main() {
  kepler_benchmarks::Nijenhuis<double> solver;
  kepler_benchmarks::run_general_benchmark("nijenhuis.csv", solver);
  return 0;
}