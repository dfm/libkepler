#ifndef KEPLER_STARTERS_HPP
#define KEPLER_STARTERS_HPP

#include <cmath>

namespace kepler {
namespace starters {

template <class T>
inline T basic_starter(const T& eccentricity, const T& mean_anomaly) noexcept {
  return mean_anomaly + 0.85 * eccentricity;
}

// https://ui.adsabs.harvard.edu/abs/1987CeMec..40..329M/abstract
template <class T>
inline T mikkola(const T& eccentricity, const T& mean_anomaly) noexcept {
  auto factor = 1. / (4. * eccentricity + 0.5);
  auto alpha = (1. - eccentricity) * factor;
  auto beta = 0.5 * mean_anomaly * factor;
  auto z = std::cbrt(beta + std::copysign(std::sqrt(beta * beta + alpha * alpha * alpha), beta));
  auto s = z - alpha / z;
  s -= 0.078 * std::pow(s, 5) / (1. + eccentricity);
  return mean_anomaly + eccentricity * s * (3. - 4. * s * s);
}

// https://ui.adsabs.harvard.edu/abs/1995CeMDA..63..101M/abstract
template <class T>
inline T markley(const T& eccentricity, const T& mean_anomaly) noexcept {
  auto m2 = mean_anomaly * mean_anomaly;
  auto ome = 1. - eccentricity;

  auto alpha = (constants::pi<T>() - mean_anomaly) / (1. + eccentricity);
  alpha *= constants::markley_factor2<T>();
  alpha += constants::markley_factor1<T>();

  auto d = 3. * ome + alpha * eccentricity;
  alpha *= d;

  auto r = mean_anomaly * (3. * alpha * (d - ome) + m2);
  auto q = 2. * alpha * ome - m2;
  auto q2 = q * q;

  auto w = std::cbrt(std::abs(r) + std::sqrt(q2 * q + r * r));
  w *= w;

  auto denom = w * (w + q) + q2;
  return (2. * r * w / denom + mean_anomaly) / d;
}

template <class T>
inline T rpp_singular(const T& eccentricity, const T& mean_anomaly) noexcept {
  auto ome = 1. - eccentricity;
  auto sqrt_ome = std::sqrt(ome);
  auto chi = mean_anomaly / (ome * sqrt_ome);
  auto lambda = std::sqrt(8. + 9. * chi * chi);
  auto s = std::cbrt(lambda + 3. * chi);
  s *= s;
  auto sigma = 6. * chi / (2. + s + 4. / s);
  auto s2 = sigma * sigma;
  auto s4 = s2 * s2;
  auto denom = 1. / (s2 + 2.);
  auto E = 1. + s2 * ome * denom *
                    ((s2 + 20.) / 60. +
                     s2 * ome * denom * denom * (s2 * s4 + 25. * s4 + 340. * s2 + 840.) / 1400.);
  return sigma * sqrt_ome * E;
}

}  // namespace starters
}  // namespace kepler

#endif