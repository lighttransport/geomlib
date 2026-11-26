/*
 * geomlib - Continuous geodesic via MMP window propagation
 *
 * Mitchell-Mount-Papadimitriou style continuous Dijkstra: windows of reachable
 * points on mesh edges are unfolded across faces and propagated in increasing
 * distance order, producing a true surface distance rather than an edge-path
 * upper bound.
 */

#ifndef GEOMLIB_GEODESIC_MMP_H
#define GEOMLIB_GEODESIC_MMP_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Shortest continuous surface distance between two vertices. Returns
 * GEOMLIB_MESH_SUCCESS and writes `length` on success. */
int geomlib_geodesic_mmp(const geomlib_trimesh *mesh,
                         uint32_t source, uint32_t target, double *length);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_GEODESIC_MMP_H */
