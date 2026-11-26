/*
 * geomlib - Ellipse point queries
 *
 * Implicit value and containment test for an oriented ellipse in 2D.
 */

#ifndef GEOMLIB_ELLIPSE_H
#define GEOMLIB_ELLIPSE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Normalized implicit value (x/a)^2 + (y/b)^2 in the ellipse frame; <= 1
 * means inside or on the boundary. */
double geomlib_ellipse_implicit(geomlib_vec2 q, geomlib_vec2 center,
                                float a, float b, float angle);

/* Returns 1 when q lies inside or on the ellipse. */
int geomlib_point_in_ellipse(geomlib_vec2 q, geomlib_vec2 center,
                             float a, float b, float angle);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ELLIPSE_H */
