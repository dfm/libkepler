#include <iostream>

#include "kepler_benchmarks/benchmark.h"
#include "kepler_benchmarks/contour.h"

int main() {
  kepler_benchmarks::Contour<8, double> solver;
  kepler_benchmarks::run_fixed_eccentricity_benchmark("contour8_fixed.csv", solver);
  return 0;
}