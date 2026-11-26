/*
 * geomlib - Approximate convex decomposition
 *
 * Clusters mesh triangles into planar regions with VSA, then computes the
 * convex hull of each region's vertices, producing a set of convex pieces.
 */

#ifndef GEOMLIB_CONVEX_DECOMPOSITION_H
#define GEOMLIB_CONVEX_DECOMPOSITION_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_trimesh *pieces;
    size_t piece_count;
} geomlib_convex_decomposition;

geomlib_convex_decomposition geomlib_convex_decomposition_create(void);
void geomlib_convex_decomposition_free(geomlib_convex_decomposition *d);

/*
 * Approximate `mesh` with `region_count` convex hull pieces.
 */
int geomlib_approx_convex_decomposition(const geomlib_trimesh *mesh,
                                        size_t region_count,
                                        geomlib_convex_decomposition *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CONVEX_DECOMPOSITION_H */
