#ifndef KEPLER_SOLVE_HPP
#define KEPLER_SOLVE_HPP

#include <cmath>
#include <cstddef>
#include <type_traits>

#include "constants.hpp"
#include "reduction.hpp"
#include "refiners.hpp"
#include "starters.hpp"

namespace kepler {

template <typename A, typename B>
struct value_type {
  static_assert(std::is_same<typename A::value_type, typename B::value_type>::value,
                "value_type of A and B must be the same");
  typedef typename A::value_type type;
};

template <typename Starter, typename Refiner>
inline void solve_one(const typename value_type<Starter, Refiner>::type& eccentricity,
                      const typename value_type<Starter, Refiner>::type& mean_anomaly,
                      typename value_type<Starter, Refiner>::type& eccentric_anomaly,
                      typename value_type<Starter, Refiner>::type& sin_eccentric_anomaly,
                      typename value_type<Starter, Refiner>::type& cos_eccentric_anomaly,
                      const Refiner& refiner, const Starter& starter) {
  using T = typename value_type<Starter, Refiner>::type;
  auto abs_mean_anom = std::abs(mean_anomaly);
  auto sgn = std::copysign(T(1.), mean_anomaly);
  T mean_anom_reduc;
  bool high = range_reduce(abs_mean_anom, mean_anom_reduc);
  auto ecc_anom_reduc = starter.start(mean_anom_reduc);
  T s, c;
  ecc_anom_reduc = refiners::refine_with_eccentricity<Refiner>::refine(
      refiner, eccentricity, mean_anom_reduc, ecc_anom_reduc, &s, &c);
  if (high) {
    eccentric_anomaly = sgn * (constants::twopi<T>() - ecc_anom_reduc);
    sin_eccentric_anomaly = -sgn * s;
    cos_eccentric_anomaly = c;
  } else {
    eccentric_anomaly = sgn * ecc_anom_reduc;
    sin_eccentric_anomaly = sgn * s;
    cos_eccentric_anomaly = c;
  }
}

template <typename Starter, typename Refiner>
inline void solve(const typename value_type<Starter, Refiner>::type& eccentricity,
                  std::size_t size,
                  const typename value_type<Starter, Refiner>::type* mean_anomaly,
                  typename value_type<Starter, Refiner>::type* eccentric_anomaly,
                  typename value_type<Starter, Refiner>::type* sin_eccentric_anomaly,
                  typename value_type<Starter, Refiner>::type* cos_eccentric_anomaly,
                  const Refiner& refiner = Refiner()) {
  const Starter starter(eccentricity);
  for (std::size_t i = 0; i < size; ++i) {
    solve_one(eccentricity, mean_anomaly[i], eccentric_anomaly[i], sin_eccentric_anomaly[i],
              cos_eccentric_anomaly[i], refiner, starter);
  }
}

template <typename Starter, typename Refiner, typename Tag = xs::unaligned_mode>
inline void solve_simd(const typename value_type<Starter, Refiner>::type& eccentricity,
                       std::size_t size,
                       const typename value_type<Starter, Refiner>::type* mean_anomaly,
                       typename value_type<Starter, Refiner>::type* eccentric_anomaly,
                       typename value_type<Starter, Refiner>::type* sin_eccentric_anomaly,
                       typename value_type<Starter, Refiner>::type* cos_eccentric_anomaly,
                       const Refiner& refiner = Refiner()) {
  using T = typename value_type<Starter, Refiner>::type;
  using B = xs::batch<T>;
  constexpr std::size_t simd_size = B::size;
  std::size_t vec_size = size - size % simd_size;
  const Starter starter(eccentricity);

  for (std::size_t i = 0; i < vec_size; i += simd_size) {
    auto mean_anom = xs::load(&(mean_anomaly[i]), Tag());
    auto sgn = xs::copysign(B(1.), mean_anom);
    auto abs_mean_anom = xs::abs(mean_anom);
    B mean_anom_reduc;
    auto high = range_reduce(abs_mean_anom, mean_anom_reduc);
    auto ecc_anom_reduc = starter.start(mean_anom_reduc);
    B s, c;
    ecc_anom_reduc = refiners::refine_with_eccentricity<Refiner>::refine(
        refiner, eccentricity, mean_anom_reduc, ecc_anom_reduc, &s, &c);
    auto ecc_anom = sgn * xs::select(high, constants::twopi<T>() - ecc_anom_reduc, ecc_anom_reduc);
    auto sin_ecc_anom = sgn * s * xs::select(high, B(-1.), B(1.));
    auto cos_ecc_anom = c;
    ecc_anom.store(&eccentric_anomaly[i], Tag());
    sin_ecc_anom.store(&sin_eccentric_anomaly[i], Tag());
    cos_ecc_anom.store(&cos_eccentric_anomaly[i], Tag());
  }

  for (std::size_t i = vec_size; i < size; ++i) {
    solve_one(eccentricity, mean_anomaly[i], eccentric_anomaly[i], sin_eccentric_anomaly[i],
              cos_eccentric_anomaly[i], refiner, starter);
  }
}

}  // namespace kepler

#endif