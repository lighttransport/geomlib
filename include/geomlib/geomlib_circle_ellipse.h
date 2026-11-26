/*
 * geomlib - Circle / ellipse overlap
 *
 * Tests overlap by projecting the circle centre onto the ellipse boundary and
 * checking centre containment.
 */

#ifndef GEOMLIB_CIRCLE_ELLIPSE_H
#define GEOMLIB_CIRCLE_ELLIPSE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int intersects;        /* 1 on overlap or containment */
    int center_inside;     /* 1 when the circle centre is inside the ellipse */
    float distance;        /* centre-to-boundary distance minus radius */
    geomlib_vec2 closest;  /* closest ellipse boundary point to the centre */
} geomlib_circle_ellipse_result;

int geomlib_circle_ellipse_overlap(geomlib_vec2 circle_center, float radius,
                                   geomlib_vec2 ellipse_center,
                                   float a, float b, float angle,
                                   geomlib_circle_ellipse_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CIRCLE_ELLIPSE_H */
