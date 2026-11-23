/*
 * geomlib - Geodesic Distance Computation
 * Based on "The Heat Method for Distance Computation" (Crane et al. 2012)
 */

#ifndef GEOMLIB_GEODESIC_H
#define GEOMLIB_GEODESIC_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Geodesic distance data */
typedef struct geomlib_geodesic_distances {
    float *distances;        /* Distance from source(s) to each vertex */
    size_t vertex_count;
    uint32_t *source_vertices;  /* Source vertex indices */
    size_t source_count;
} geomlib_geodesic_distances;

/* Create geodesic distance structure */
geomlib_geodesic_distances geomlib_geodesic_create(size_t vertex_count);

/* Free geodesic distance structure */
void geomlib_geodesic_free(geomlib_geodesic_distances *geodesic);

/* Heat-method geodesic distances (Crane et al. 2012): solves the heat equation,
 * normalizes its gradient, and solves a Poisson equation. A time_step <= 0 selects
 * mean_edge_length^2. Falls back to edge-following Dijkstra only if the linear
 * solves fail. */
int geomlib_compute_geodesic_heat(
    const geomlib_trimesh *mesh,
    const uint32_t *source_vertices,
    size_t source_count,
    float time_step,
    geomlib_geodesic_distances *distances
);

/* Compute geodesic distance from a single source vertex */
int geomlib_compute_geodesic_from_vertex(
    const geomlib_trimesh *mesh,
    uint32_t source_vertex,
    geomlib_geodesic_distances *distances
);

/* Compute all-pairs geodesic distances (expensive!) */
int geomlib_compute_geodesic_all_pairs(
    const geomlib_trimesh *mesh,
    float **distance_matrix  /* Output: vertex_count x vertex_count matrix */
);

/* Helper: Compute mean edge length for mesh */
float geomlib_compute_mean_edge_length_for_geodesic(const geomlib_trimesh *mesh);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_GEODESIC_H */
