#ifndef KEPLER_SOLVE_HPP
#define KEPLER_SOLVE_HPP

#include <cmath>
#include <cstddef>
#include <type_traits>

#include "./constants.hpp"
#include "./reduction.hpp"
#include "./refiners.hpp"
#include "./starters.hpp"

namespace kepler {

template <typename A, typename B>
struct value_type {
  static_assert(std::is_same<typename A::value_type, typename B::value_type>::value,
                "value_type of A and B must be the same");
  typedef typename A::value_type type;
};

template <typename Starter, typename Refiner>
inline void solve(const typename value_type<Starter, Refiner>::type& eccentricity,
                  std::size_t size,
                  const typename value_type<Starter, Refiner>::type* mean_anomaly,
                  typename value_type<Starter, Refiner>::type* eccentric_anomaly,
                  const Refiner& refiner = Refiner()) {
  using T = typename value_type<Starter, Refiner>::type;
  const Starter starter(eccentricity);
  for (std::size_t i = 0; i < size; ++i) {
    auto mean_anom_reduc = std::abs(mean_anomaly[i]);
    bool high = range_reduce(mean_anom_reduc, mean_anom_reduc);
    auto ecc_anom_reduc = starter.start(mean_anom_reduc);
    ecc_anom_reduc = refiner.refine(eccentricity, mean_anom_reduc, ecc_anom_reduc);
    eccentric_anomaly[i] = std::copysign(
        high ? constants::twopi<T>() - ecc_anom_reduc : ecc_anom_reduc, mean_anomaly[i]);
  }
}

template <typename Starter, typename Refiner, typename Tag = xs::unaligned_mode>
inline void solve_simd(const typename value_type<Starter, Refiner>::type& eccentricity,
                       std::size_t size,
                       const typename value_type<Starter, Refiner>::type* mean_anomaly,
                       typename value_type<Starter, Refiner>::type* eccentric_anomaly,
                       const Refiner& refiner = Refiner()) {
  using T = typename value_type<Starter, Refiner>::type;
  using B = xs::batch<T>;
  constexpr std::size_t simd_size = B::size;
  std::size_t vec_size = size - size % simd_size;
  const Starter starter(eccentricity);

  for (std::size_t i = 0; i < vec_size; i += simd_size) {
    auto mean_anom = xs::load(&(mean_anomaly[i]), Tag());
    auto abs_mean_anom = xs::abs(mean_anom);
    B mean_anom_reduc;
    auto high = range_reduce(abs_mean_anom, mean_anom_reduc);
    auto ecc_anom_reduc = starter.start(mean_anom_reduc);
    ecc_anom_reduc = refiner.refine(eccentricity, mean_anom_reduc, ecc_anom_reduc);
    auto ecc_anom = xs::copysign(
        xs::select(high, constants::twopi<T>() - ecc_anom_reduc, ecc_anom_reduc), mean_anom);
    ecc_anom.store(&eccentric_anomaly[i], Tag());
  }

  for (std::size_t i = vec_size; i < size; ++i) {
    auto mean_anom_reduc = std::abs(mean_anomaly[i]);
    bool high = range_reduce(mean_anom_reduc, mean_anom_reduc);
    auto ecc_anom_reduc = starter.start(mean_anom_reduc);
    ecc_anom_reduc = refiner.refine(eccentricity, mean_anom_reduc, ecc_anom_reduc);
    eccentric_anomaly[i] = std::copysign(
        high ? constants::twopi<T>() - ecc_anom_reduc : ecc_anom_reduc, mean_anomaly[i]);
  }
}

}  // namespace kepler

#endif