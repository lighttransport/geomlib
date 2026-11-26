/*
 * geomlib - Point-cloud local curvature estimation
 *
 * Computes a local covariance from k nearest neighbours and reports the
 * normalized smallest eigenvalue as a curvature measure, with its eigenvector
 * as the local normal.
 */

#ifndef GEOMLIB_POINT_CURVATURE_H
#define GEOMLIB_POINT_CURVATURE_H

#include "geomlib_surface.h"

#ifdef __cplusplus
extern "C" {
#endif

/* `curvature` and `normals` are optional output arrays of cloud->point_count
 * entries. `k` <= 0 selects 16 neighbours. */
int geomlib_pointcloud_curvature(const geomlib_pointcloud *cloud, int k,
                                 float *curvature, geomlib_vec3 *normals);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POINT_CURVATURE_H */
