#ifndef KEPLER_SOLVE_HPP
#define KEPLER_SOLVE_HPP

#include <cstddef>

#include "./refiners.hpp"
#include "./starters.hpp"

namespace kepler {

template <typename Refiner, typename Starter = typename Refiner::default_starter,
          typename T = typename Refiner::value_type>
inline void solve(const T& eccentricity, std::size_t size, const T* mean_anomaly,
                  T* eccentric_anomaly, const Refiner& refiner = Refiner()) {
  const Starter starter(eccentricity);
  for (std::size_t i = 0; i < size; ++i) {
    auto initial = starter.start(mean_anomaly[i]);
    eccentric_anomaly[i] = refiner.refine(eccentricity, mean_anomaly[i], initial);
  }
}

template <typename T>
inline void solve_non_iterative_1(const T& eccentricity, std::size_t size, const T* mean_anomaly,
                                  T* eccentric_anomaly) {
  solve<refiners::non_iterative<4, T>, starters::markley<T>>(eccentricity, size, mean_anomaly,
                                                             eccentric_anomaly);
}

}  // namespace kepler

#endif