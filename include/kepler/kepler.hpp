// Copyright 2021-2024 Dan Foreman-Mackey
//
// Distributed under the terms of the Apache 2.0 License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef KEPLER_KEPLER_HPP
#define KEPLER_KEPLER_HPP

#include <cstdint>

#include "kepler/kepler/refiners.hpp"
#include "kepler/kepler/solver.hpp"
#include "kepler/kepler/starters.hpp"

namespace kepler {

template <typename T>
void solve(std::size_t size, const T* eccentricity, std::size_t batch_size, const T* mean_anomaly,
           T* eccentric_anomaly, T* sin_eccentric_anomaly, T* cos_eccentric_anomaly) {
  for (std::size_t n = 0; n < size; ++n) {
    solver::solve_simd<starters::raposo_pulido_brandt<T>, refiners::brandt<T>>(
        eccentricity[n], batch_size, mean_anomaly, eccentric_anomaly, sin_eccentric_anomaly,
        cos_eccentric_anomaly);
    mean_anomaly += batch_size;
    eccentric_anomaly += batch_size;
    sin_eccentric_anomaly += batch_size;
    cos_eccentric_anomaly += batch_size;
  }
}

}  // namespace kepler
#endif