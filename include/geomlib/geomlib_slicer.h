/*
 * geomlib - Mesh slicing
 *
 * Intersects a triangle mesh with a plane and returns the resulting polylines.
 */

#ifndef GEOMLIB_SLICER_H
#define GEOMLIB_SLICER_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 *points;      /* Concatenated polyline points */
    size_t point_count;
    uint32_t *starts;          /* Start index of each polyline */
    uint32_t *counts;          /* Number of points in each polyline */
    size_t polyline_count;
} geomlib_slice_result;

geomlib_slice_result geomlib_slice_result_create(void);
void geomlib_slice_result_free(geomlib_slice_result *r);

/* Slice a mesh by the plane dot(normal, x - point) = 0. */
int geomlib_slice_mesh(const geomlib_trimesh *mesh,
                       geomlib_vec3 point,
                       geomlib_vec3 normal,
                       float weld_eps,
                       geomlib_slice_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SLICER_H */
