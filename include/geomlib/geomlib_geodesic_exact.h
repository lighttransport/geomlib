/*
 * geomlib - Continuous geodesic shortest-path length via unfolding
 *
 * Unfolds the triangle strip along an edge-graph shortest path into a plane,
 * so the straight-line distance between the unfolded endpoints is the exact
 * shortest path within that face sequence (the continuous Dijkstra/MMP
 * local-optimal solution for the strip).
 */

#ifndef GEOMLIB_GEODESIC_EXACT_H
#define GEOMLIB_GEODESIC_EXACT_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Exact length of the shortest surface path from `source` to `target`,
 * confined to the face sequence of an edge-graph shortest path. */
int geomlib_geodesic_length_unfolding(const geomlib_trimesh *mesh,
                                      uint32_t source,
                                      uint32_t target,
                                      double *length);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_GEODESIC_EXACT_H */
