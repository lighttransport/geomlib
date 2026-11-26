/*
 * geomlib - Exact graph-geodesic diameter of a triangle mesh
 *
 * Finds the pair of vertices with the largest shortest-path distance.
 */

#ifndef GEOMLIB_GEODESIC_DIAMETER_H
#define GEOMLIB_GEODESIC_DIAMETER_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t source, target;
    float distance;
} geomlib_geodesic_diameter;

/* Exhaustive vertex-pair search via Dijkstra from every vertex. */
int geomlib_geodesic_diameter_exact(const geomlib_trimesh *mesh,
                                    geomlib_geodesic_diameter *out);

/* Faster approximate search: start from a vertex, then from its farthest
 * vertex (double sweep). */
int geomlib_geodesic_diameter_sweep(const geomlib_trimesh *mesh,
                                    geomlib_geodesic_diameter *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_GEODESIC_DIAMETER_H */
