#ifndef _KEPLER_BENCHMARKS_ITERATIVE_H_
#define _KEPLER_BENCHMARKS_ITERATIVE_H_

#include "math_utils.h"
#include "solver.h"

namespace kepler_benchmarks {

// A reference first order iterative solver extracted from "batman":
// https://github.com/lkreidberg/batman/blob/master/c_src/_rsky.c
//
// From that source: calculates the eccentric anomaly (see Seager Exoplanets
// book:  Murray & Correia eqn. 5 -- see section 3)
template <typename T>
struct FirstOrderRef : public Solver<T> {
  int max_iter = 30;
  T tolerance = 1e-7;

  INLINE_OR_DEVICE T compute_eccentric_anomaly(const T& mean_anomaly,
                                               const T& eccentricity) override {
    T E = mean_anomaly + sign(sin(mean_anomaly)) * 0.85 * eccentricity;
    T fe, fs;
    for (int i = 0; i < max_iter; ++i) {
      fe = E - eccentricity * sin(E) - mean_anomaly;
      if (fabs(fe) < tolerance) break;
      fs = 1 - eccentricity * cos(E);
      E = E - fe / fs;
    }
    return E;
  }
};

// A reference third order iterative solver extracted from "radvel":
// https://github.com/California-Planet-Search/radvel/blob/master/src/kepler.c
template <typename T>
struct ThirdOrderRef : public Solver<T> {
  int max_iter = 30;
  T tolerance = 1e-12;

  INLINE_OR_DEVICE T compute_eccentric_anomaly(const T& mean_anomaly,
                                               const T& eccentricity) override {
    T E = mean_anomaly + sign(sin(mean_anomaly)) * 0.85 * eccentricity;
    T fi, fip, fipp, fippp, d1;

    for (int i = 0; i < max_iter; ++i) {
      fi = (E - eccentricity * sin(E) - mean_anomaly);
      if (fabs(fi) < tolerance) break;

      fip = 1 - eccentricity * cos(E);
      fipp = eccentricity * sin(E);
      fippp = 1 - fip;

      d1 = -fi / fip;
      d1 = -fi / (fip + d1 * fipp / 2.0);
      d1 = -fi / (fip + d1 * fipp / 2.0 + d1 * d1 * fippp / 6.0);
      E += d1;
    }

    return E;
  }
};

}  // namespace kepler_benchmarks

#endif
