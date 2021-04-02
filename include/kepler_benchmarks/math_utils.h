#ifndef _KEPLER_BENCHMARKS_MATH_UTILS_H_
#define _KEPLER_BENCHMARKS_MATH_UTILS_H_

#include <cmath>

template <typename T>
inline int sign(T x) {
  return (x > 0) - (x < 0);
}

#endif