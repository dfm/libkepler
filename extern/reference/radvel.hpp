#ifndef KEPLER_EXTERN_RADVEL_HPP
#define KEPLER_EXTERN_RADVEL_HPP

#include <cmath>

namespace kepler {
namespace reference {

template <typename T>
inline int sign(T x) {
  return (x > 0) - (x < 0);
}

// A reference third order iterative solver extracted from "radvel":
// https://github.com/California-Planet-Search/radvel/blob/master/src/kepler.c
struct radvel {
  typedef double value_type;
  int max_iterations = 30;
  double tolerance = 1e-12;
  double eccentricity = 0.0;

  radvel() {}
  radvel(double tolerance) : tolerance(tolerance) {}
  radvel(int max_iterations, double tolerance)
      : max_iterations(max_iterations), tolerance(tolerance) {}

  inline void setup(const double& eccentricity) { this->eccentricity = eccentricity; }

  inline double solve(const double& mean_anomaly) const {
    double E = mean_anomaly + sign(std::sin(mean_anomaly)) * 0.85 * eccentricity;
    double fi, fip, fipp, fippp, d1;

    for (int i = 0; i < max_iterations; ++i) {
      fi = (E - eccentricity * sin(E) - mean_anomaly);
      if (std::abs(fi) < tolerance) break;

      fip = 1 - eccentricity * std::cos(E);
      fipp = eccentricity * std::sin(E);
      fippp = 1 - fip;

      d1 = -fi / fip;
      d1 = -fi / (fip + d1 * fipp / 2.0);
      d1 = -fi / (fip + d1 * fipp / 2.0 + d1 * d1 * fippp / 6.0);
      E += d1;
    }

    return E;
  }
};

}  // namespace reference
}  // namespace kepler

#endif
