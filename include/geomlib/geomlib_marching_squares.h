/*
 * geomlib - 2D marching squares contour extraction
 *
 * Extracts zero-level contours from a scalar grid as polylines, with
 * consistent saddle resolution and shared grid-edge vertices.
 */

#ifndef GEOMLIB_MARCHING_SQUARES_H
#define GEOMLIB_MARCHING_SQUARES_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 *points;
    size_t point_count;
    uint32_t *starts;
    uint32_t *counts;
    size_t polyline_count;
} geomlib_contour_result;

geomlib_contour_result geomlib_contour_result_create(void);
void geomlib_contour_result_free(geomlib_contour_result *r);

/* Extract zero contours from a row-major nx x ny scalar grid with spacing h
 * and origin. Values are sampled at grid nodes (not cell centers). */
int geomlib_marching_squares(const float *values, size_t nx, size_t ny,
                             float h, geomlib_vec2 origin,
                             geomlib_contour_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MARCHING_SQUARES_H */
