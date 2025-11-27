/*
 * geomlib - Mesh Decimation using Quadric Error Metrics
 * Based on "Surface Simplification Using Quadric Error Metrics" (Garland & Heckbert, SIGGRAPH 1997)
 */

#ifndef GEOMLIB_DECIMATION_H
#define GEOMLIB_DECIMATION_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Quadric error matrix (4x4 symmetric matrix, stored as 10 values) */
typedef struct geomlib_quadric {
    double a2, ab, ac, ad;  /* Row 0: a^2, ab, ac, ad */
    double b2, bc, bd;      /* Row 1: b^2, bc, bd */
    double c2, cd;          /* Row 2: c^2, cd */
    double d2;              /* Row 3: d^2 */
} geomlib_quadric;

/* Decimation parameters */
typedef struct geomlib_decimation_params {
    float target_reduction;     /* 0.0 to 1.0, fraction of vertices to remove */
    size_t target_vertex_count; /* Alternative: target number of vertices (0 = use reduction) */
    size_t target_triangle_count; /* Alternative: target number of triangles (0 = use reduction/vertices) */

    float max_error;            /* Maximum allowed quadric error (0 = no limit) */
    bool preserve_boundaries;   /* Preserve mesh boundaries */
    bool preserve_topology;     /* Prevent non-manifold geometry */
    float boundary_weight;      /* Weight multiplier for boundary edges (default: 1000.0) */
} geomlib_decimation_params;

/* Decimation statistics */
typedef struct geomlib_decimation_stats {
    size_t original_vertices;
    size_t original_triangles;
    size_t final_vertices;
    size_t final_triangles;
    size_t edges_collapsed;
    float reduction_ratio;
    float max_error;
    float mean_error;
} geomlib_decimation_stats;

/* Create default decimation parameters */
geomlib_decimation_params geomlib_decimation_params_default(float target_reduction);

/* Decimate mesh using QEM
 * Returns GEOMLIB_MESH_SUCCESS on success
 * The input mesh is modified in place
 */
int geomlib_decimate_qem(
    geomlib_trimesh *mesh,
    const geomlib_decimation_params *params,
    geomlib_decimation_stats *stats  /* Optional, can be NULL */
);

/* Simplified decimation interface - reduce to target percentage */
int geomlib_decimate_to_ratio(
    geomlib_trimesh *mesh,
    float target_ratio  /* 0.0 to 1.0, fraction of vertices to keep */
);

/* Simplified decimation interface - reduce to target vertex count */
int geomlib_decimate_to_vertex_count(
    geomlib_trimesh *mesh,
    size_t target_vertices
);

/* Simplified decimation interface - reduce to target triangle count */
int geomlib_decimate_to_triangle_count(
    geomlib_trimesh *mesh,
    size_t target_triangles
);

/* Helper functions */

/* Initialize a quadric from a plane (ax + by + cz + d = 0) */
geomlib_quadric geomlib_quadric_from_plane(float a, float b, float c, float d);

/* Add two quadrics */
geomlib_quadric geomlib_quadric_add(geomlib_quadric q1, geomlib_quadric q2);

/* Compute error for a vertex position */
double geomlib_quadric_error(geomlib_quadric q, geomlib_vec3 v);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_DECIMATION_H */
