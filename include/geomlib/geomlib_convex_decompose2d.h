/*
 * geomlib - 2D polygon convex decomposition
 *
 * Recursively cuts a simple polygon along valid diagonals from reflex
 * vertices, choosing the cut that minimizes new reflex vertices.
 */

#ifndef GEOMLIB_CONVEX_DECOMPOSE2D_H
#define GEOMLIB_CONVEX_DECOMPOSE2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 *parts;     /* concatenated polygons */
    size_t *offsets;         /* part_count+1 offsets into parts */
    size_t part_count;
} geomlib_poly_decomposition;

geomlib_poly_decomposition geomlib_poly_decomposition_create(void);
void geomlib_poly_decomposition_free(geomlib_poly_decomposition *d);

int geomlib_convex_decompose2d(const geomlib_vec2 *polygon, size_t n,
                               geomlib_poly_decomposition *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CONVEX_DECOMPOSE2D_H */
