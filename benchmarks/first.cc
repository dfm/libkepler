#include <iostream>

#include "kepler_benchmarks/benchmark.h"
#include "kepler_benchmarks/iterative.h"

int main() {
  kepler_benchmarks::FirstOrderRef<double> first;
  kepler_benchmarks::run_general_benchmark("first.csv", first);
  return 0;
}