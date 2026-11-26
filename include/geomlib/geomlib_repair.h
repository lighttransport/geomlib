/*
 * geomlib - Robust topology repair
 *
 * Clean-room implementations of conservative, topology-only repairs used to
 * turn imperfect input into a manifold triangle mesh.
 */

#ifndef GEOMLIB_REPAIR_H
#define GEOMLIB_REPAIR_H

#include "geomlib_mesh.h"

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Split edges with more than two incident faces by detaching one of the
 * incident triangles at a time (all three of its vertices are duplicated).
 * This is conservative: it is guaranteed to reduce edge over-use and never
 * changes any triangle's geometry. Returns the number of detached triangles.
 */
int geomlib_repair_nonmanifold_edges(geomlib_trimesh *mesh,
                                     size_t *detached_out);

/*
 * Build a triangle mesh from a polygon soup. Points closer than `weld_eps`
 * are merged, degenerate and duplicate triangles are removed, and winding
 * is made consistent across shared edges.
 */
int geomlib_soup_to_mesh(const geomlib_vec3 *points, size_t point_count,
                         const uint32_t *indices, size_t index_count,
                         float weld_eps, geomlib_trimesh *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_REPAIR_H */
