/*
 * geomlib - Heat kernel signature (HKS) shape descriptor
 *
 * HKS_i(t) = sum_k exp(-lambda_k t) phi_k(i)^2, computed from the lowest
 * Laplace-Beltrami eigenpairs.
 */

#ifndef GEOMLIB_HKS_H
#define GEOMLIB_HKS_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Compute HKS for each vertex at each time. `out` is row-major with shape
 * (vertex_count, ntimes); `times` must be positive and ascending. */
int geomlib_compute_hks(const geomlib_trimesh *mesh, size_t eigen_count,
                        const double *times, size_t ntimes, double *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_HKS_H */
