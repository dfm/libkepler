#ifndef KEPLER_STARTERS_HPP
#define KEPLER_STARTERS_HPP

#include <cmath>

#include "constants.hpp"
#include "simd.hpp"

namespace kepler {
namespace starters {

template <typename T>
struct noop {
  typedef T value_type;
  noop(T) {}

  inline T start(const T& mean_anomaly) const { return mean_anomaly; }

  template <typename A>
  inline xs::batch<T, A> start(const xs::batch<T, A>& mean_anomaly) const {
    return mean_anomaly;
  }
};

template <typename T>
struct basic {
  typedef T value_type;
  T eccentricity;
  basic(T eccentricity) : eccentricity(eccentricity) {}

  inline T start(const T& mean_anomaly) const { return mean_anomaly + T(0.85) * eccentricity; }

  template <typename A>
  inline xs::batch<T, A> start(const xs::batch<T, A>& mean_anomaly) const {
    return mean_anomaly + xs::batch<T, A>(T(0.85) * eccentricity);
  }
};

// https://ui.adsabs.harvard.edu/abs/1987CeMec..40..329M/abstract
template <typename T>
struct mikkola {
  typedef T value_type;
  T eccentricity, factor, alpha, alpha3;
  mikkola(T eccentricity)
      : eccentricity(eccentricity),
        factor(T(1.) / (T(4.) * eccentricity + T(0.5))),
        alpha((T(1.) - eccentricity) * factor),
        alpha3(alpha * alpha * alpha) {}

  inline T start(const T& mean_anomaly) const {
    auto beta = T(0.5) * mean_anomaly * factor;
    auto z = std::cbrt(beta + std::copysign(std::sqrt(beta * beta + alpha3), beta));
    auto s = z - alpha / z;
    s -= T(0.078) * std::pow(s, T(5.)) / (T(1.) + eccentricity);
    return mean_anomaly + eccentricity * s * (T(3.) - T(4.) * s * s);
  }

  template <typename A>
  inline xs::batch<T, A> start(const xs::batch<T, A>& mean_anomaly) const {
    using B = xs::batch<T, A>;
    auto beta = B(T(0.5) * factor) * mean_anomaly;
    auto z = xs::cbrt(beta + xs::copysign(xs::sqrt(xs::fma(beta, beta, B(alpha3))), beta));
    auto s = xs::fnma(B(alpha), T(1.) / z, z);
    s -= B(T(0.078) / (T(1.) + eccentricity)) * xs::pow(s, 5);
    return xs::fma(B(eccentricity) * s, xs::fnma(B(4.) * s, s, B(3.)), mean_anomaly);
  }
};

// https://ui.adsabs.harvard.edu/abs/1995CeMDA..63..101M/abstract
template <typename T>
struct markley {
  typedef T value_type;
  T eccentricity;
  markley(T eccentricity) : eccentricity(eccentricity) {}

  inline T start(const T& mean_anomaly) const {
    auto m2 = mean_anomaly * mean_anomaly;
    auto ome = T(1.) - eccentricity;

    auto alpha = (constants::pi<T>() - mean_anomaly) / (T(1.) + eccentricity);
    alpha *= constants::markley_factor2<T>();
    alpha += constants::markley_factor1<T>();

    auto d = T(3.) * ome + alpha * eccentricity;
    alpha *= d;

    auto r = mean_anomaly * (T(3.) * alpha * (d - ome) + m2);
    auto q = T(2.) * alpha * ome - m2;
    auto q2 = q * q;

    auto w = std::cbrt(std::abs(r) + std::sqrt(q2 * q + r * r));
    w *= w;

    auto denom = w * (w + q) + q2;
    return (T(2.) * r * w / denom + mean_anomaly) / d;
  }

  template <typename A>
  inline xs::batch<T, A> start(const xs::batch<T, A>& mean_anomaly) const {
    using B = xs::batch<T, A>;
    auto m2 = mean_anomaly * mean_anomaly;
    auto ome = T(1.) - eccentricity;

    auto alpha = xs::fma(B(constants::markley_factor2<T>() / (T(1.) + eccentricity)),
                         constants::pi<T>() - mean_anomaly, B(constants::markley_factor1<T>()));

    auto d = xs::fma(B(eccentricity), alpha, B(T(3.) * ome));
    alpha *= d;

    auto r = mean_anomaly * xs::fma(B(T(3.)) * alpha, d - B(ome), m2);
    auto q = xs::fms(B(T(2.) * ome), alpha, m2);
    auto q2 = q * q;

    auto w = xs::cbrt(xs::abs(r) + xs::sqrt(xs::fma(q2, q, r * r)));
    w *= w;

    auto denom = xs::fma(w, w + q, q2);
    return xs::fma(B(T(2.)) * r / denom, w, mean_anomaly) / d;
  }
};

// https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.1702R/abstract
// https://ui.adsabs.harvard.edu/abs/2021AJ....162..186B/abstract
template <typename T>
struct raposo_pulido_brandt {
  typedef T value_type;
  T eccentricity, ome, sqrt_ome, bounds[13], table[78];

  raposo_pulido_brandt(T eccentricity)
      : eccentricity(eccentricity), ome(T(1.) - eccentricity), sqrt_ome(std::sqrt(ome)) {
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
    table[1] = T(1.) / (T(1.) - eccentricity);
    table[2] = T(0.);

    x = T(1.) / (T(1.) - g2c_e);
    table[7] = x;
    table[8] = -T(0.5) * g2s_e * x * x * x;

    x = T(1.) / (T(1.) - g3c_e);
    table[13] = x;
    table[14] = -T(0.5) * g3s_e * x * x * x;

    x = T(1.) / (T(1.) - g4c_e);
    table[19] = x;
    table[20] = -T(0.5) * g4s_e * x * x * x;

    x = T(1.) / (T(1.) - g5c_e);
    table[25] = x;
    table[26] = -T(0.5) * g5s_e * x * x * x;

    x = T(1.) / (T(1.) - g6c_e);
    table[31] = x;
    table[32] = -T(0.5) * g6s_e * x * x * x;

    table[37] = T(1.);
    table[38] = -T(0.5) * eccentricity;

    x = T(1.) / (T(1.) + g6c_e);
    table[43] = x;
    table[44] = -T(0.5) * g6s_e * x * x * x;

    x = T(1.) / (T(1.) + g5c_e);
    table[49] = x;
    table[50] = -T(0.5) * g5s_e * x * x * x;

    x = T(1.) / (T(1.) + g4c_e);
    table[55] = x;
    table[56] = -T(0.5) * g4s_e * x * x * x;

    x = T(1.) / (T(1.) + g3c_e);
    table[61] = x;
    table[62] = -T(0.5) * g3s_e * x * x * x;

    x = T(1.) / (T(1.) + g2c_e);
    table[67] = x;
    table[68] = -T(0.5) * g2s_e * x * x * x;

    table[73] = T(1.) / (T(1.) + eccentricity);
    table[74] = T(0.);

    for (int i = 0; i < 12; i++) {
      int k = 6 * i;
      table[k] = T(i) * constants::pio12<T>();

      auto idx = T(1.) / (bounds[i + 1] - bounds[i]);
      auto B0 = idx * (-table[k + 2] - idx * (table[k + 1] - idx * constants::pio12<T>()));
      auto B1 = idx * (-T(2.) * table[k + 2] - idx * (table[k + 1] - table[k + 7]));
      auto B2 = idx * (table[k + 8] - table[k + 2]);

      table[k + 3] = B2 - T(4.) * B1 + T(10.) * B0;
      table[k + 4] = (-T(2.) * B2 + T(7.) * B1 - T(15.) * B0) * idx;
      table[k + 5] = (B2 - T(3.) * B1 + T(6.) * B0) * idx * idx;
    }
  }

  inline T singular(const T& mean_anomaly) const {
    auto chi = mean_anomaly / (ome * sqrt_ome);
    auto lambda = std::sqrt(T(8.) + T(9.) * chi * chi);
    auto s = std::cbrt(lambda + T(3.) * chi);
    s *= s;
    auto sigma = T(6.) * chi / (T(2.) + s + T(4.) / s);
    auto s2 = sigma * sigma;
    auto denom = T(1.) / (s2 + T(2.));
    auto arg = s2 * ome * denom * denom * (s2 * (s2 * (s2 + T(25.)) + T(340.)) + T(840.));
    auto E = T(1.) + s2 * ome * denom * ((s2 + T(20.)) / T(60.) + arg / T(1400.));
    return sigma * sqrt_ome * E;
  }

  template <typename A>
  inline xs::batch<T, A> singular(const xs::batch<T, A>& mean_anomaly) const {
    using B = xs::batch<T, A>;
    auto chi = mean_anomaly / (ome * sqrt_ome);
    auto lambda = xs::sqrt(xs::fma(B(T(9.)) * chi, chi, B(T(8.))));
    auto s = xs::cbrt(xs::fma(B(T(3.)), chi, lambda));
    s *= s;
    auto sigma = B(T(6.)) * chi / (B(T(2.)) + s + B(T(4.)) / s);
    auto s2 = sigma * sigma;
    auto denom = B(T(1.)) / (s2 + B(T(2.)));
    auto arg = B(ome) * s2 * denom * denom *
               xs::fma(s2, xs::fma(s2, (s2 + T(25.)), B(T(340.))), B(T(840.)));
    auto E = xs::fma(B(ome) * s2, denom * ((s2 + T(20.)) / T(60.) + arg / T(1400.)), B(T(1.)));
    return sigma * sqrt_ome * E;
  }

  inline T lookup(const T& mean_anomaly) const {
    int j;
    for (j = 11; j > 0; --j)
      if (mean_anomaly > bounds[j]) break;
    auto k = 6 * j;
    auto dx = mean_anomaly - bounds[j];
    return math::horner(dx, table[k], table[k + 1], table[k + 2], table[k + 3], table[k + 4],
                        table[k + 5]);
  }

  template <typename A>
  inline xs::batch<T, A> lookup(const xs::batch<T, A>& mean_anomaly) const {
    using B = xs::batch<T, A>;
    using I = typename xs::as_integer_t<B>;
    static_assert(B::size == I::size, "integer batch size must match float batch size");

    // We fall back on the serial algorithm for the lookup table because it
    // seems to be slower to use a vectorized version.
    std::array<typename I::value_type, I::size> idx;
    std::array<typename B::value_type, B::size> val;
    mean_anomaly.store_aligned(val.data());
    for (size_t n = 0; n < I::size; ++n) {
      auto v = val[n];
      typename I::value_type j = 0;
      for (j = 11; j > 0; --j)
        if (v > bounds[j]) break;
      idx[n] = j;
    }

    // Now back to vectorized code
    auto j = xs::load_aligned(idx.data());
    auto k = I(6) * j;
    auto dx = mean_anomaly - B::gather(bounds, j);
    return math::horner(dx, B::gather(table, k), B::gather(table, k + 1), B::gather(table, k + 2),
                        B::gather(table, k + 3), B::gather(table, k + 4), B::gather(table, k + 5));
  }

  inline T start(const T& mean_anomaly) const {
    if ((eccentricity < T(0.78)) || (T(2.) * mean_anomaly + ome > T(0.2))) {
      return lookup(mean_anomaly);
    } else {
      return singular(mean_anomaly);
    }
  }

  template <typename A>
  inline xs::batch<T, A> start(const xs::batch<T, A>& mean_anomaly) const {
    using B = xs::batch<T, A>;
    auto flag =
        (B(eccentricity) < B(T(0.78))) | (xs::fma(B(T(2.)), mean_anomaly, B(ome)) > B(T(0.2)));
    auto fastpath = lookup(mean_anomaly);
    if (xs::all(flag)) return fastpath;
    return xs::select(flag, fastpath, singular(mean_anomaly));
  }
};

}  // namespace starters
}  // namespace kepler

#endif