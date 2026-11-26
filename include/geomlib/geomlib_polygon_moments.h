/*
 * geomlib - 2D polygon area moments
 *
 * Computes signed area, centroid and second moments about the centroid.
 */

#ifndef GEOMLIB_POLYGON_MOMENTS_H
#define GEOMLIB_POLYGON_MOMENTS_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double area;      /* absolute area */
    geomlib_vec2 centroid;
    double ixx, ixy, iyy;   /* second moments about the centroid */
} geomlib_polygon_moments;

int geomlib_compute_polygon_moments(const geomlib_vec2 *polygon, size_t n,
                            geomlib_polygon_moments *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYGON_MOMENTS_H */
