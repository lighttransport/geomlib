/*
 * geomlib - Exact graph geodesic distances (Dijkstra)
 *
 * Computes shortest-path distances along mesh edges from one or more source
 * vertices. Uses `geomlib_geodesic_distances` for the result.
 */

#ifndef GEOMLIB_GEODESIC_DIJKSTRA_H
#define GEOMLIB_GEODESIC_DIJKSTRA_H

#include "geomlib_geodesic.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Multi-source Dijkstra over the triangle-edge graph. */
int geomlib_geodesic_dijkstra(const geomlib_trimesh *mesh,
                              const uint32_t *sources, size_t source_count,
                              geomlib_geodesic_distances *out);

/* Shorter helper for a single source. */
int geomlib_geodesic_dijkstra_from_vertex(const geomlib_trimesh *mesh,
                                          uint32_t source,
                                          geomlib_geodesic_distances *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_GEODESIC_DIJKSTRA_H */
