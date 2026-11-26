/*
 * geomlib - 2D circle / line and circle / segment intersection
 *
 * Analytic quadratic along a parametric line or segment.
 */

#ifndef GEOMLIB_CIRCLE_LINE_H
#define GEOMLIB_CIRCLE_LINE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    GEOMLIB_CIRCLE_LINE_NONE = 0,
    GEOMLIB_CIRCLE_LINE_TANGENT,
    GEOMLIB_CIRCLE_LINE_TWO_POINTS
} geomlib_circle_line_kind;

typedef struct {
    geomlib_circle_line_kind kind;
    geomlib_vec2 p0, p1;   /* intersection points (equal when tangent) */
} geomlib_circle_line_intersection;

/* Intersection of a circle with the infinite line through a and b. */
geomlib_circle_line_intersection geomlib_circle_line_intersect(
    geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 center, float radius);

/* Intersection of a circle with the finite segment [a,b]. */
geomlib_circle_line_intersection geomlib_circle_segment_intersect(
    geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 center, float radius);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CIRCLE_LINE_H */
