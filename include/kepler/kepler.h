// Copyright 2021-2024 Dan Foreman-Mackey
//
// Distributed under the terms of the Apache 2.0 License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef KEPLER_KEPLER_H
#define KEPLER_KEPLER_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void kepler_solve(size_t size, const double* eccentricity, size_t batch_size,
                  const double* mean_anomaly, double* eccentric_anomaly,
                  double* sin_eccentric_anomaly, double* cos_eccentric_anomaly);
void kepler_solvef(size_t size, const float* eccentricity, size_t batch_size,
                   const float* mean_anomaly, float* eccentric_anomaly,
                   float* sin_eccentric_anomaly, float* cos_eccentric_anomaly);

#ifdef __cplusplus
}
#endif
#endif