/*
 * geomlib - Remeshing Module
 *
 * Implements the Botsch-Kobbelt incremental remeshing algorithm for
 * creating uniform or adaptive triangle meshes with high-quality elements.
 *
 * References:
 * - "A Remeshing Approach to Multiresolution Modeling"
 *   Botsch & Kobbelt, SGP 2004
 * - PMP Library: https://www.pmp-library.org/remeshing.html
 * - Geometry Central: https://geometry-central.net/surface/algorithms/remeshing/
 *
 * Algorithm:
 * The algorithm iteratively applies local mesh modifications:
 * 1. Edge Split: Split edges longer than 4/3 * target_length
 * 2. Edge Collapse: Collapse edges shorter than 4/5 * target_length
 * 3. Edge Flip: Flip edges to improve vertex valence (target = 6)
 * 4. Tangential Smoothing: Smooth vertices along tangent plane
 * 5. Project to surface: Optional projection to preserve shape
 */

#ifndef GEOMLIB_REMESHING_H
#define GEOMLIB_REMESHING_H

#include "geomlib_mesh.h"
#include "geomlib_linalg.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Remeshing mode */
typedef enum geomlib_remeshing_mode {
    GEOMLIB_REMESH_UNIFORM,      /* Uniform target edge length */
    GEOMLIB_REMESH_ADAPTIVE      /* Adaptive based on curvature */
} geomlib_remeshing_mode;

/* Remeshing parameters */
typedef struct geomlib_remeshing_params {
    /* Mode */
    geomlib_remeshing_mode mode;

    /* Target edge length (for uniform mode) */
    float target_edge_length;

    /* Min/max edge length (for adaptive mode) */
    float min_edge_length;
    float max_edge_length;

    /* Number of iterations */
    int iterations;

    /* Feature preservation */
    float feature_angle;         /* Dihedral angle threshold in degrees */
    bool preserve_boundaries;    /* Preserve mesh boundaries */

    /* Projection */
    bool project_to_surface;     /* Project vertices back to original surface */

    /* Advanced parameters */
    float split_threshold;       /* Edge split if length > split_threshold * target */
    float collapse_threshold;    /* Edge collapse if length < collapse_threshold * target */
    int smoothing_iterations;    /* Tangential smoothing iterations per remesh iteration */
} geomlib_remeshing_params;

/* Default parameters */
geomlib_remeshing_params geomlib_remeshing_params_default(
    geomlib_remeshing_mode mode,
    float target_edge_length
);

/* Main remeshing function */
int geomlib_remesh(
    geomlib_trimesh *mesh,
    const geomlib_remeshing_params *params
);

/* Uniform remeshing (simplified interface) */
int geomlib_remesh_uniform(
    geomlib_trimesh *mesh,
    float target_edge_length,
    int iterations
);

/* Adaptive remeshing (simplified interface) */
int geomlib_remesh_adaptive(
    geomlib_trimesh *mesh,
    float min_edge_length,
    float max_edge_length,
    int iterations
);

/* Compute target edge length from desired vertex count */
float geomlib_compute_target_edge_length(
    const geomlib_trimesh *mesh,
    size_t target_vertex_count
);

/* Compute average edge length of mesh */
float geomlib_compute_mean_edge_length(const geomlib_trimesh *mesh);

/* Statistics */
typedef struct geomlib_remeshing_stats {
    int edges_split;
    int edges_collapsed;
    int edges_flipped;
    float min_edge_length;
    float max_edge_length;
    float mean_edge_length;
    float min_angle;            /* Minimum triangle angle in degrees */
    float max_angle;            /* Maximum triangle angle in degrees */
    int irregular_vertices;     /* Vertices with valence != 6 */
} geomlib_remeshing_stats;

/* Compute mesh quality statistics */
void geomlib_compute_remeshing_stats(
    const geomlib_trimesh *mesh,
    geomlib_remeshing_stats *stats
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_REMESHING_H */
