#ifndef KEPLER_HOUSEHOLDER_HPP
#define KEPLER_HOUSEHOLDER_HPP

#include <xsimd/xsimd.hpp>

#include "./constants.hpp"
#include "./math.hpp"

namespace xs = xsimd;

namespace kepler {
namespace detail {

template <typename T>
inline T householder_eval(const T& f0, const T& f1) {
  return -f0 / f1;
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2) {
  auto d = householder_eval(f0, f1);
  return -f0 / (f1 + d * constants::hh2<T>() * f2);
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2, const T& f3) {
  auto d = householder_eval(f0, f1, f2);
  return -f0 / (f1 + d * (constants::hh2<T>() * f2 + d * constants::hh3<T>() * f3));
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2, const T& f3, const T& f4) {
  auto d = householder_eval(f0, f1, f2, f3);
  return -f0 / (f1 + d * (constants::hh2<T>() * f2 +
                          d * (constants::hh3<T>() * f3 + d * constants::hh4<T>() * f4)));
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2, const T& f3, const T& f4,
                          const T& f5) {
  auto d = householder_eval(f0, f1, f2, f3);
  return -f0 / (f1 + d * (constants::hh2<T>() * f2 +
                          d * (constants::hh3<T>() * f3 +
                               d * (constants::hh4<T>() * f4 + d * constants::hh5<T>() * f5))));
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2, const T& f3, const T& f4,
                          const T& f5, const T& f6) {
  auto d = householder_eval(f0, f1, f2, f3);
  return -f0 /
         (f1 + d * (constants::hh2<T>() * f2 +
                    d * (constants::hh3<T>() * f3 +
                         d * (constants::hh4<T>() * f4 +
                              d * (constants::hh5<T>() * f5 + d * constants::hh6<T>() * f6)))));
}

template <typename T>
inline T householder_eval(const T& f0, const T& f1, const T& f2, const T& f3, const T& f4,
                          const T& f5, const T& f6, const T& f7) {
  auto d = householder_eval(f0, f1, f2, f3);
  return -f0 /
         (f1 +
          d * (constants::hh2<T>() * f2 +
               d * (constants::hh3<T>() * f3 +
                    d * (constants::hh4<T>() * f4 +
                         d * (constants::hh5<T>() * f5 +
                              d * (constants::hh6<T>() * f6 + d * constants::hh7<T>() * f7))))));
}

template <int order>
struct householder {
  template <typename T>
  static inline T step(const T& eccentricity, const T& mean_anomaly, T& eccentric_anomaly);

  template <typename A, typename T>
  static inline xs::batch<T, A> step(const T& eccentricity, const xs::batch<T, A>& mean_anomaly,
                                     xs::batch<T, A>& eccentric_anomaly);
};

template <>
struct householder<1> {
  template <typename T>
  static inline T step(const T& eccentricity, const T& mean_anomaly, T& eccentric_anomaly) {
    T s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = T(1.) - eccentricity;
    auto f0 = eccentricity * s + eccentric_anomaly * ome - mean_anomaly;
    auto f1 = eccentricity * c + ome;
    eccentric_anomaly += householder_eval(f0, f1);
    return f0;
  }

  template <typename A, typename T>
  static inline xs::batch<T, A> step(const T& eccentricity, const xs::batch<T, A>& mean_anomaly,
                                     xs::batch<T, A>& eccentric_anomaly) {
    using B = xs::batch<T, A>;
    B s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = B(T(1.) - eccentricity);
    auto f0 = xs::fma(B(eccentricity), s, eccentric_anomaly * ome) - mean_anomaly;
    auto f1 = xs::fma(B(eccentricity), c, B(ome));
    eccentric_anomaly += householder_eval(f0, f1);
    return f0;
  }
};

template <>
struct householder<2> {
  template <typename T>
  static inline T step(const T& eccentricity, const T& mean_anomaly, T& eccentric_anomaly) {
    T s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = T(1.) - eccentricity;
    auto f0 = eccentricity * s + eccentric_anomaly * ome - mean_anomaly;
    auto f1 = eccentricity * c + ome;
    auto f2 = eccentricity * (eccentric_anomaly - s);
    eccentric_anomaly += householder_eval(f0, f1, f2);
    return f0;
  }

  template <typename A, typename T>
  static inline xs::batch<T, A> step(const T& eccentricity, const xs::batch<T, A>& mean_anomaly,
                                     xs::batch<T, A>& eccentric_anomaly) {
    using B = xs::batch<T, A>;
    B s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = B(T(1.) - eccentricity);
    auto f0 = xs::fma(B(eccentricity), s, eccentric_anomaly * ome) - mean_anomaly;
    auto f1 = xs::fma(B(eccentricity), c, B(ome));
    auto f2 = B(eccentricity) * (eccentric_anomaly - s);
    eccentric_anomaly += householder_eval(f0, f1, f2);
    return f0;
  }
};

template <>
struct householder<3> {
  template <typename T>
  static inline T step(const T& eccentricity, const T& mean_anomaly, T& eccentric_anomaly) {
    T s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = T(1.) - eccentricity;
    auto f0 = eccentricity * s + eccentric_anomaly * ome - mean_anomaly;
    auto f1 = eccentricity * c + ome;
    auto f2 = eccentricity * (eccentric_anomaly - s);
    auto f3 = T(1.) - f1;
    eccentric_anomaly += householder_eval(f0, f1, f2, f3);
    return f0;
  }

  template <typename A, typename T>
  static inline xs::batch<T, A> step(const T& eccentricity, const xs::batch<T, A>& mean_anomaly,
                                     xs::batch<T, A>& eccentric_anomaly) {
    using B = xs::batch<T, A>;
    B s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = B(T(1.) - eccentricity);
    auto f0 = xs::fma(B(eccentricity), s, eccentric_anomaly * ome) - mean_anomaly;
    auto f1 = xs::fma(B(eccentricity), c, B(ome));
    auto f2 = B(eccentricity) * (eccentric_anomaly - s);
    auto f3 = B(1.) - f1;
    eccentric_anomaly += householder_eval(f0, f1, f2, f3);
    return f0;
  }
};

template <>
struct householder<4> {
  template <typename T>
  static inline T step(const T& eccentricity, const T& mean_anomaly, T& eccentric_anomaly) {
    T s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = T(1.) - eccentricity;
    auto f0 = eccentricity * s + eccentric_anomaly * ome - mean_anomaly;
    auto f1 = eccentricity * c + ome;
    auto f2 = eccentricity * (eccentric_anomaly - s);
    auto f3 = T(1.) - f1;
    auto f4 = -f2;
    eccentric_anomaly += householder_eval(f0, f1, f2, f3, f4);
    return f0;
  }

  template <typename A, typename T>
  static inline xs::batch<T, A> step(const T& eccentricity, const xs::batch<T, A>& mean_anomaly,
                                     xs::batch<T, A>& eccentric_anomaly) {
    using B = xs::batch<T, A>;
    B s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = B(T(1.) - eccentricity);
    auto f0 = xs::fma(B(eccentricity), s, eccentric_anomaly * ome) - mean_anomaly;
    auto f1 = xs::fma(B(eccentricity), c, B(ome));
    auto f2 = B(eccentricity) * (eccentric_anomaly - s);
    auto f3 = B(1.) - f1;
    auto f4 = -f2;
    eccentric_anomaly += householder_eval(f0, f1, f2, f3, f4);
    return f0;
  }
};

template <>
struct householder<5> {
  template <typename T>
  static inline T step(const T& eccentricity, const T& mean_anomaly, T& eccentric_anomaly) {
    T s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = T(1.) - eccentricity;
    auto f0 = eccentricity * s + eccentric_anomaly * ome - mean_anomaly;
    auto f1 = eccentricity * c + ome;
    auto f2 = eccentricity * (eccentric_anomaly - s);
    auto f3 = T(1.) - f1;
    auto f4 = -f2;
    auto f5 = -f3;
    eccentric_anomaly += householder_eval(f0, f1, f2, f3, f4, f5);
    return f0;
  }

  template <typename A, typename T>
  static inline xs::batch<T, A> step(const T& eccentricity, const xs::batch<T, A>& mean_anomaly,
                                     xs::batch<T, A>& eccentric_anomaly) {
    using B = xs::batch<T, A>;
    B s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = B(T(1.) - eccentricity);
    auto f0 = xs::fma(B(eccentricity), s, eccentric_anomaly * ome) - mean_anomaly;
    auto f1 = xs::fma(B(eccentricity), c, B(ome));
    auto f2 = B(eccentricity) * (eccentric_anomaly - s);
    auto f3 = B(1.) - f1;
    auto f4 = -f2;
    auto f5 = -f3;
    eccentric_anomaly += householder_eval(f0, f1, f2, f3, f4, f5);
    return f0;
  }
};

template <>
struct householder<6> {
  template <typename T>
  static inline T step(const T& eccentricity, const T& mean_anomaly, T& eccentric_anomaly) {
    T s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = T(1.) - eccentricity;
    auto f0 = eccentricity * s + eccentric_anomaly * ome - mean_anomaly;
    auto f1 = eccentricity * c + ome;
    auto f2 = eccentricity * (eccentric_anomaly - s);
    auto f3 = T(1.) - f1;
    auto f4 = -f2;
    auto f5 = -f3;
    auto f6 = -f4;
    eccentric_anomaly += householder_eval(f0, f1, f2, f3, f4, f5, f6);
    return f0;
  }

  template <typename A, typename T>
  static inline xs::batch<T, A> step(const T& eccentricity, const xs::batch<T, A>& mean_anomaly,
                                     xs::batch<T, A>& eccentric_anomaly) {
    using B = xs::batch<T, A>;
    B s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = B(T(1.) - eccentricity);
    auto f0 = xs::fma(B(eccentricity), s, eccentric_anomaly * ome) - mean_anomaly;
    auto f1 = xs::fma(B(eccentricity), c, B(ome));
    auto f2 = B(eccentricity) * (eccentric_anomaly - s);
    auto f3 = B(1.) - f1;
    auto f4 = -f2;
    auto f5 = -f3;
    auto f6 = -f4;
    eccentric_anomaly += householder_eval(f0, f1, f2, f3, f4, f5, f6);
    return f0;
  }
};

template <>
struct householder<7> {
  template <typename T>
  static inline T step(const T& eccentricity, const T& mean_anomaly, T& eccentric_anomaly) {
    T s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = T(1.) - eccentricity;
    auto f0 = eccentricity * s + eccentric_anomaly * ome - mean_anomaly;
    auto f1 = eccentricity * c + ome;
    auto f2 = eccentricity * (eccentric_anomaly - s);
    auto f3 = T(1.) - f1;
    auto f4 = -f2;
    auto f5 = -f3;
    auto f6 = -f4;
    auto f7 = -f5;
    eccentric_anomaly += householder_eval(f0, f1, f2, f3, f4, f5, f6, f7);
    return f0;
  }

  template <typename A, typename T>
  static inline xs::batch<T, A> step(const T& eccentricity, const xs::batch<T, A>& mean_anomaly,
                                     xs::batch<T, A>& eccentric_anomaly) {
    using B = xs::batch<T, A>;
    B s, c;
    ::kepler::detail::sin_cos_reduc(eccentric_anomaly, s, c);
    auto ome = B(T(1.) - eccentricity);
    auto f0 = xs::fma(B(eccentricity), s, eccentric_anomaly * ome) - mean_anomaly;
    auto f1 = xs::fma(B(eccentricity), c, B(ome));
    auto f2 = B(eccentricity) * (eccentric_anomaly - s);
    auto f3 = B(1.) - f1;
    auto f4 = -f2;
    auto f5 = -f3;
    auto f6 = -f4;
    auto f7 = -f5;
    eccentric_anomaly += householder_eval(f0, f1, f2, f3, f4, f5, f6, f7);
    return f0;
  }
};

}  // namespace detail
}  // namespace kepler

#endif