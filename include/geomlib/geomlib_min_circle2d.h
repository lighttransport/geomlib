/*
 * geomlib - Exact 2D minimum enclosing circle (Welzl's algorithm)
 */

#ifndef GEOMLIB_MIN_CIRCLE2D_H
#define GEOMLIB_MIN_CIRCLE2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 center;
    float radius;
} geomlib_min_circle2d;

int geomlib_minimum_enclosing_circle(const geomlib_vec2 *points, size_t count,
                                     geomlib_min_circle2d *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MIN_CIRCLE2D_H */
