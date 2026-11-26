/*
 * geomlib - Minimum-area bounding rectangle (rotating calipers)
 *
 * Finds the enclosing rectangle of smallest area around a planar point set.
 */

#ifndef GEOMLIB_MIN_AREA_RECT_H
#define GEOMLIB_MIN_AREA_RECT_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 center;
    geomlib_vec2 axes[2];      /* unit edge directions */
    geomlib_vec2 extents;     /* half-widths along the axes */
    float area;
} geomlib_min_area_rect;

int geomlib_minimum_area_rectangle(const geomlib_vec2 *points, size_t n,
                                   geomlib_min_area_rect *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MIN_AREA_RECT_H */
