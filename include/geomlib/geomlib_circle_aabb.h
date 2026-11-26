/*
 * geomlib - 2D circle / axis-aligned rectangle overlap
 *
 * Clamps the circle centre to the rectangle and compares with the radius.
 */

#ifndef GEOMLIB_CIRCLE_AABB_H
#define GEOMLIB_CIRCLE_AABB_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;         /* 1 on overlap or tangency */
    float distance;         /* centre-to-rectangle distance minus radius */
    geomlib_vec2 closest;   /* closest point on the rectangle */
} geomlib_circle_aabb_result;

int geomlib_circle_aabb_overlap(geomlib_vec2 center, float radius,
                                geomlib_vec2 bmin, geomlib_vec2 bmax,
                                geomlib_circle_aabb_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CIRCLE_AABB_H */
