/*
 * geomlib - 2D line / line intersection
 *
 * Infinite-line intersection with parallel and collinear classification.
 */

#ifndef GEOMLIB_LINE2D_H
#define GEOMLIB_LINE2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    GEOMLIB_LINE2D_NONE = 0,
    GEOMLIB_LINE2D_POINT,
    GEOMLIB_LINE2D_COLLINEAR
} geomlib_line2d_kind;

typedef struct {
    geomlib_line2d_kind kind;
    geomlib_vec2 point;   /* unique intersection point */
    double t;             /* parameter along the first line */
    double u;             /* parameter along the second line */
} geomlib_line2d_intersection;

geomlib_line2d_intersection geomlib_line2d_intersect(
    geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c, geomlib_vec2 d);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_LINE2D_H */
