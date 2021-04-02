#include <iostream>

#include "kepler_benchmarks/benchmark.h"
#include "kepler_benchmarks/iterative.h"

int main() {
  kepler_benchmarks::ThirdOrderRef<double> third;
  kepler_benchmarks::run_general_benchmark("third.csv", third);
  return 0;
}