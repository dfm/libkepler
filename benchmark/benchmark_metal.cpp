#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <vector>

#include "kepler/metal/metal.hpp"

TEST_CASE("Metal", "[metal]") {
  NS::AutoreleasePool *pool = NS::AutoreleasePool::alloc()->init();

  MTL::Device *device = MTL::CreateSystemDefaultDevice();
  if (!device) {
    std::cerr << "Failed to find device.\n";
    exit(-1);
  }

  using T = float;

  const std::size_t size = 10000;
  std::vector<T> mean_anomaly(size), ecc_anomaly(size);
  for (size_t m = 0; m < size; ++m) {
    mean_anomaly[m] = T(100.) * m / T(size - 1) - T(50.);
  }

  kepler::metal::solver<T> solver(
      device, "/Users/dforemanmackey/src/dfm/kepler.cpp/build/src/kepler/metal/default.metallib");

  BENCHMARK("Set data") { return solver.set_data(0.5, size, mean_anomaly.data()); };

  solver.set_data(0.5, size, mean_anomaly.data());
  BENCHMARK("Execute") { return solver.execute(); };

  solver.execute();
  BENCHMARK("Set data") { return solver.get_results(ecc_anomaly.data()); };

  device->release();
  pool->release();
}