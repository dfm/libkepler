#ifndef KEPLER_EXTERN_BATMAN_HPP
#define KEPLER_EXTERN_BATMAN_HPP

#include <cmath>

namespace kepler {
namespace reference {

// A reference first order iterative solver extracted from "batman":
// https://github.com/lkreidberg/batman/blob/master/c_src/_rsky.c
//
// From that source: calculates the eccentric anomaly (see Seager Exoplanets
// book:  Murray & Correia eqn. 5 -- see section 3)
struct batman {
  typedef double value_type;
  int max_iterations = 30;
  double tolerance = 1e-7;
  double eccentricity = 0.0;

  batman() {}
  batman(double tolerance) : tolerance(tolerance) {}
  batman(int max_iterations, double tolerance)
      : max_iterations(max_iterations), tolerance(tolerance) {}

  inline void setup(const double& eccentricity) { this->eccentricity = eccentricity; }

  inline double solve(const double& mean_anomaly) const {
    double E = mean_anomaly;
    double fe, fs;
    for (int i = 0; i < max_iterations; ++i) {
      fe = E - eccentricity * std::sin(E) - mean_anomaly;
      if (std::abs(fe) < tolerance) break;
      fs = 1 - eccentricity * std::cos(E);
      E = E - fe / fs;
    }
    return E;
  }
};

}  // namespace reference
}  // namespace kepler

#endif
