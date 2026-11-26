/*
 * geomlib - Mesh Parameterization Module
 *
 * Implements mesh UV parameterization algorithms for mapping 3D surfaces to 2D.
 *
 * References:
 * - "Least Squares Conformal Maps for Automatic Texture Atlas Generation"
 *   Lévy, Petitjean, Ray, Maillot, SIGGRAPH 2002
 * - "A Local/Global Approach to Mesh Parameterization"
 *   Liu et al. (ARAP parameterization)
 * - "Parameterization of Faceted Surfaces for Meshing using Angle-Based Flattening"
 *   Sheffer & de Sturler (ABF/ABF++)
 * - CGAL Surface Mesh Parameterization documentation
 *
 * Algorithms Implemented:
 * - LSCM: Least Squares Conformal Maps (conformal/angle-preserving)
 * - Tutte/Harmonic: Convex combination mapping (guaranteed bijective)
 */

#ifndef GEOMLIB_PARAMETERIZATION_H
#define GEOMLIB_PARAMETERIZATION_H

#include "geomlib_mesh.h"
#include "geomlib_linalg.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Parameterization method */
typedef enum geomlib_param_method {
    GEOMLIB_PARAM_LSCM,        /* Least Squares Conformal Maps (angle-preserving) */
    GEOMLIB_PARAM_HARMONIC,    /* Harmonic/Tutte (guaranteed bijective) */
    GEOMLIB_PARAM_ARAP         /* As-Rigid-As-Possible (shape-preserving, future) */
} geomlib_param_method;

/* Boundary shape for parameterization */
typedef enum geomlib_boundary_shape {
    GEOMLIB_BOUNDARY_CIRCLE,   /* Map boundary to unit circle */
    GEOMLIB_BOUNDARY_SQUARE,   /* Map boundary to unit square */
    GEOMLIB_BOUNDARY_NATURAL   /* Use natural boundary (LSCM free boundary) */
} geomlib_boundary_shape;

/* Parameterization parameters */
typedef struct geomlib_param_params {
    geomlib_param_method method;
    geomlib_boundary_shape boundary;

    /* For methods requiring boundary pinning */
    uint32_t pin_vertex_0;     /* First pinned vertex (use UINT32_MAX for auto) */
    uint32_t pin_vertex_1;     /* Second pinned vertex (use UINT32_MAX for auto) */

    /* Solver parameters */
    int max_iterations;        /* Max iterations for iterative solvers */
    float tolerance;           /* Convergence tolerance */

    /* Normalization */
    bool normalize_uvs;        /* Normalize UVs to [0,1] range */
} geomlib_param_params;

/* UV coordinates storage */
typedef struct geomlib_uv_coords {
    geomlib_vec2 *uvs;         /* UV coordinates per vertex */
    size_t vertex_count;
    float *distortion;         /* Optional: per-triangle distortion metric */
    size_t triangle_count;
} geomlib_uv_coords;

/* Default parameters */
geomlib_param_params geomlib_param_params_default(geomlib_param_method method);

/* Main parameterization function */
int geomlib_parameterize(
    const geomlib_trimesh *mesh,
    const geomlib_param_params *params,
    geomlib_uv_coords *uvs
);

/* LSCM parameterization (conformal mapping) */
int geomlib_parameterize_lscm(
    const geomlib_trimesh *mesh,
    geomlib_boundary_shape boundary,
    geomlib_uv_coords *uvs
);

/* Harmonic/Tutte parameterization (guaranteed bijective) */
int geomlib_parameterize_harmonic(
    const geomlib_trimesh *mesh,
    geomlib_boundary_shape boundary,
    geomlib_uv_coords *uvs
);

/* ARAP parameterization (Liu et al. 2008 local/global solver).
 * Boundary vertices are pinned to a circle (or square). */
int geomlib_parameterize_arap(
    const geomlib_trimesh *mesh,
    geomlib_boundary_shape boundary,
    int max_iterations,
    geomlib_uv_coords *uvs
);

/* LSCM with a free boundary: only two vertices are pinned (at (0,0) and
 * (1,0)), all other boundary vertices are solved for. */
int geomlib_parameterize_lscm_free(
    const geomlib_trimesh *mesh,
    uint32_t pin_vertex_0,
    uint32_t pin_vertex_1,
    geomlib_uv_coords *uvs
);

/* Helper: Find boundary loop of mesh */
int geomlib_find_boundary_loop(
    const geomlib_trimesh *mesh,
    uint32_t **boundary_vertices,
    size_t *boundary_count
);

/* Helper: Pin boundary vertices to circle/square */
void geomlib_pin_boundary_circle(
    const uint32_t *boundary_vertices,
    size_t boundary_count,
    geomlib_vec2 *uvs
);

void geomlib_pin_boundary_square(
    const uint32_t *boundary_vertices,
    size_t boundary_count,
    geomlib_vec2 *uvs
);

/* Compute distortion metrics */
void geomlib_compute_param_distortion(
    const geomlib_trimesh *mesh,
    const geomlib_uv_coords *uvs
);

/* Normalize UVs to [0,1] range */
void geomlib_normalize_uvs(geomlib_uv_coords *uvs);

/* UV coordinates management */
geomlib_uv_coords geomlib_uv_coords_create(size_t vertex_count, size_t triangle_count);
void geomlib_uv_coords_free(geomlib_uv_coords *uvs);

/* Statistics */
typedef struct geomlib_param_stats {
    float angle_distortion;    /* RMS angle distortion */
    float area_distortion;     /* RMS area distortion */
    float stretch_L2;          /* L2 stretch metric */
    float stretch_Linf;        /* L∞ stretch metric */
    bool has_flips;            /* Whether any triangles are flipped */
    int flipped_count;         /* Number of flipped triangles */
} geomlib_param_stats;

void geomlib_compute_param_stats(
    const geomlib_trimesh *mesh,
    const geomlib_uv_coords *uvs,
    geomlib_param_stats *stats
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_PARAMETERIZATION_H */
