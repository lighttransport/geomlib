/*
 * geomlib - Largest inscribed circle of a simple polygon
 *
 * Approximates the Chebyshev centre by evaluating Voronoi vertices of
 * densely sampled boundary points and keeping the one furthest from the
 * polygon edges.
 */

#ifndef GEOMLIB_INSCRIBED_CIRCLE_H
#define GEOMLIB_INSCRIBED_CIRCLE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 center;
    float radius;
} geomlib_inscribed_circle;

/* Compute the largest inscribed circle. `spacing` <= 0 selects 0.02 * bbox
 * diagonal. Returns GEOMLIB_MESH_SUCCESS on success. */
int geomlib_largest_inscribed_circle(const geomlib_vec2 *polygon, size_t n,
                                     float spacing,
                                     geomlib_inscribed_circle *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_INSCRIBED_CIRCLE_H */
