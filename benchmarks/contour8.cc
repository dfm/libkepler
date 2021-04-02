#include <iostream>

#include "kepler_benchmarks/benchmark.h"
#include "kepler_benchmarks/contour.h"

int main() {
  kepler_benchmarks::Contour<8, double> solver;
  kepler_benchmarks::run_general_benchmark("contour8.csv", solver);
  return 0;
}