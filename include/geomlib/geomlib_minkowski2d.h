/*
 * geomlib - 2D convex Minkowski sum
 *
 * Minkowski sum of two convex polygons, built as the convex hull of all
 * pairwise vertex sums.
 */

#ifndef GEOMLIB_MINKOWSKI2D_H
#define GEOMLIB_MINKOWSKI2D_H

#include "geomlib_hull2d.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Compute A + B. For non-convex inputs the result is the Minkowski sum of
 * their convex hulls. Returns GEOMLIB_MESH_SUCCESS on success. */
int geomlib_minkowski_sum2d(const geomlib_vec2 *a, size_t na,
                            const geomlib_vec2 *b, size_t nb,
                            geomlib_hull2d *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MINKOWSKI2D_H */
