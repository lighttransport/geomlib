/*
 * geomlib - Minimum enclosing sphere (Welzl's algorithm)
 *
 * Computes the smallest sphere containing a point set, exactly for the
 * supported dimensions.
 */

#ifndef GEOMLIB_MINIBALL_H
#define GEOMLIB_MINIBALL_H

#include "geomlib/geomlib_bounding_sphere.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Exact smallest enclosing sphere of `count` 3D points. */
int geomlib_minimum_enclosing_sphere(const geomlib_vec3 *points, size_t count,
                                     geomlib_bounding_sphere *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MINIBALL_H */
