/*
 * geomlib - Least-squares 3D plane fitting
 *
 * Fits a plane to 3D points by PCA: the normal is the smallest eigenvector
 * of the centred covariance matrix.
 */

#ifndef GEOMLIB_PLANE_FIT_H
#define GEOMLIB_PLANE_FIT_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 point;     /* a point on the plane (the centroid) */
    geomlib_vec3 normal;    /* unit normal */
    float rms;
} geomlib_plane_fit;

int geomlib_fit_plane(const geomlib_vec3 *points, size_t count,
                      geomlib_plane_fit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_PLANE_FIT_H */
