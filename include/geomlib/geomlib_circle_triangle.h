/*
 * geomlib - 2D circle / triangle intersection
 *
 * Classifies circle-triangle overlap using containment and boundary distance.
 */

#ifndef GEOMLIB_CIRCLE_TRIANGLE_H
#define GEOMLIB_CIRCLE_TRIANGLE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;         /* 1 when the circle meets the triangle */
    int center_inside;      /* 1 when the circle centre is inside the triangle */
    int contains_triangle;  /* 1 when all three vertices lie in the circle */
    float signed_distance;  /* centre-to-boundary distance minus radius */
    geomlib_vec2 closest;   /* closest triangle boundary point to the centre */
} geomlib_circle_triangle_result;

int geomlib_circle_triangle_intersect(geomlib_vec2 center, float radius,
                                      geomlib_vec2 a, geomlib_vec2 b,
                                      geomlib_vec2 c,
                                      geomlib_circle_triangle_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CIRCLE_TRIANGLE_H */
