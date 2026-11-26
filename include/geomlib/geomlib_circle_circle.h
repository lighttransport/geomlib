/*
 * geomlib - 2D circle / circle intersection
 *
 * Analytic radical-line solution with tangency, containment and disjoint
 * classification.
 */

#ifndef GEOMLIB_CIRCLE_CIRCLE_H
#define GEOMLIB_CIRCLE_CIRCLE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    GEOMLIB_CIRCLE2D_NONE = 0,
    GEOMLIB_CIRCLE2D_TANGENT,
    GEOMLIB_CIRCLE2D_TWO_POINTS,
    GEOMLIB_CIRCLE2D_CONTAINED,
    GEOMLIB_CIRCLE2D_COINCIDENT
} geomlib_circle2d_kind;

typedef struct {
    geomlib_circle2d_kind kind;
    geomlib_vec2 p0, p1;   /* intersection points (equal when tangent) */
} geomlib_circle2d_intersection;

geomlib_circle2d_intersection geomlib_circle2d_intersect(
    geomlib_vec2 c0, float r0, geomlib_vec2 c1, float r1);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CIRCLE_CIRCLE_H */
