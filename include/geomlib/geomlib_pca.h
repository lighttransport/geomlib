/*
 * geomlib - Principal component analysis of 3D points
 *
 * Computes the mean, covariance eigenvalues (descending) and orthonormal
 * eigenvectors of a point set.
 */

#ifndef GEOMLIB_PCA_H
#define GEOMLIB_PCA_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 mean;
    double eigenvalues[3];      /* descending */
    geomlib_vec3 eigenvectors[3];
} geomlib_pca_result;

int geomlib_compute_pca(const geomlib_vec3 *points, size_t count,
                        geomlib_pca_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_PCA_H */
