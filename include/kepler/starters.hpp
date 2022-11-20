#ifndef KEPLER_STARTERS_HPP
#define KEPLER_STARTERS_HPP

#include <cmath>

#include "./constants.hpp"

namespace kepler {
namespace starters {

template <typename T>
struct noop {
  noop(T) {}
  inline T start(const T& mean_anomaly) const { return mean_anomaly; }
};

template <typename T>
struct basic {
  T eccentricity;
  basic(T eccentricity) : eccentricity(eccentricity) {}
  inline T start(const T& mean_anomaly) const { return mean_anomaly + 0.85 * eccentricity; }
};

// https://ui.adsabs.harvard.edu/abs/1987CeMec..40..329M/abstract
template <typename T>
struct mikkola {
  T eccentricity;
  mikkola(T eccentricity) : eccentricity(eccentricity) {}
  inline T start(const T& mean_anomaly) const {
    auto factor = 1. / (4. * eccentricity + 0.5);
    auto alpha = (1. - eccentricity) * factor;
    auto beta = 0.5 * mean_anomaly * factor;
    auto z = std::cbrt(beta + std::copysign(std::sqrt(beta * beta + alpha * alpha * alpha), beta));
    auto s = z - alpha / z;
    s -= 0.078 * std::pow(s, 5) / (1. + eccentricity);
    return mean_anomaly + eccentricity * s * (3. - 4. * s * s);
  }
};

// https://ui.adsabs.harvard.edu/abs/1995CeMDA..63..101M/abstract
template <typename T>
struct markley {
  T eccentricity;
  markley(T eccentricity) : eccentricity(eccentricity) {}
  inline T start(const T& mean_anomaly) const {
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
};

// https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.1702R/abstract
// https://ui.adsabs.harvard.edu/abs/2021AJ....162..186B/abstract
template <typename T>
struct rppb {
  T eccentricity, ome, sqrt_ome, bounds[13], table[78];

  rppb(T eccentricity)
      : eccentricity(eccentricity), ome(1. - eccentricity), sqrt_ome(std::sqrt(ome)) {
    auto g2s_e = constants::rppb_g2s<T>() * eccentricity;
    auto g3s_e = constants::rppb_g3s<T>() * eccentricity;
    auto g4s_e = constants::rppb_g4s<T>() * eccentricity;
    auto g5s_e = constants::rppb_g5s<T>() * eccentricity;
    auto g6s_e = constants::rppb_g6s<T>() * eccentricity;
    auto g2c_e = g6s_e;
    auto g3c_e = g5s_e;
    auto g4c_e = g4s_e;
    auto g5c_e = g3s_e;
    auto g6c_e = g2s_e;

    bounds[0] = T(0.);
    bounds[1] = constants::pio12<T>() - g2s_e;
    bounds[2] = constants::pio6<T>() - g3s_e;
    bounds[3] = constants::pio4<T>() - g4s_e;
    bounds[4] = constants::pio3<T>() - g5s_e;
    bounds[5] = constants::fivepio12<T>() - g6s_e;
    bounds[6] = constants::pio2<T>() - eccentricity;
    bounds[7] = constants::sevenpio12<T>() - g6s_e;
    bounds[8] = constants::twopio3<T>() - g5s_e;
    bounds[9] = constants::threepio4<T>() - g4s_e;
    bounds[10] = constants::fivepio6<T>() - g3s_e;
    bounds[11] = constants::elevenpio12<T>() - g2s_e;
    bounds[12] = constants::pi<T>();

    T x;
    table[1] = 1. / (1. - eccentricity);
    table[2] = T(0.);

    x = 1. / (1 - g2c_e);
    table[7] = x;
    table[8] = -0.5 * g2s_e * x * x * x;

    x = 1. / (1. - g3c_e);
    table[13] = x;
    table[14] = -0.5 * g3s_e * x * x * x;

    x = 1. / (1. - g4c_e);
    table[19] = x;
    table[20] = -0.5 * g4s_e * x * x * x;

    x = 1. / (1. - g5c_e);
    table[25] = x;
    table[26] = -0.5 * g5s_e * x * x * x;

    x = 1. / (1. - g6c_e);
    table[31] = x;
    table[32] = -0.5 * g6s_e * x * x * x;

    table[37] = T(1.);
    table[38] = -0.5 * eccentricity;

    x = 1. / (1. + g6c_e);
    table[43] = x;
    table[44] = -0.5 * g6s_e * x * x * x;

    x = 1. / (1. + g5c_e);
    table[49] = x;
    table[50] = -0.5 * g5s_e * x * x * x;

    x = 1. / (1. + g4c_e);
    table[55] = x;
    table[56] = -0.5 * g4s_e * x * x * x;

    x = 1. / (1. + g3c_e);
    table[61] = x;
    table[62] = -0.5 * g3s_e * x * x * x;

    x = 1. / (1. + g2c_e);
    table[67] = x;
    table[68] = -0.5 * g2s_e * x * x * x;

    table[73] = 1. / (1 + eccentricity);
    table[74] = T(0.);

    for (int i = 0; i < 12; i++) {
      int k = 6 * i;
      table[k] = T(i) * constants::pio12<T>();

      auto idx = 1. / (bounds[i + 1] - bounds[i]);
      auto B0 = idx * (-table[k + 2] - idx * (table[k + 1] - idx * constants::pio12<T>()));
      auto B1 = idx * (-2. * table[k + 2] - idx * (table[k + 1] - table[k + 7]));
      auto B2 = idx * (table[k + 8] - table[k + 2]);

      table[k + 3] = B2 - 4. * B1 + 10. * B0;
      table[k + 4] = (-2. * B2 + 7. * B1 - 15. * B0) * idx;
      table[k + 5] = (B2 - 3. * B1 + 6. * B0) * idx * idx;
    }
  }

  inline T singular(const T& mean_anomaly) const {
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

  inline T start(const T& mean_anomaly) const {
    if ((eccentricity < 0.78) || (2. * mean_anomaly + ome > 0.2)) {
      int j;
      for (j = 11; j > 0; --j)
        if (mean_anomaly > bounds[j]) break;
      auto k = 6 * j;
      auto dx = mean_anomaly - bounds[j];
      return table[k] +
             dx * (table[k + 1] +
                   dx * (table[k + 2] +
                         dx * (table[k + 3] + dx * (table[k + 4] + dx * table[k + 5]))));
    } else {
      return singular(mean_anomaly);
    }
  }
};

}  // namespace starters
}  // namespace kepler

#endif