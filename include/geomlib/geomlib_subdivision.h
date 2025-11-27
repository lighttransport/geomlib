/*
 * geomlib - Subdivision Surfaces
 * Loop and Catmull-Clark subdivision schemes
 */

#ifndef GEOMLIB_SUBDIVISION_H
#define GEOMLIB_SUBDIVISION_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Subdivision scheme */
typedef enum geomlib_subdivision_scheme {
    GEOMLIB_SUBDIV_LOOP,          /* Loop subdivision (triangles) */
    GEOMLIB_SUBDIV_BUTTERFLY,     /* Butterfly subdivision (interpolating) */
    GEOMLIB_SUBDIV_MIDPOINT       /* Simple midpoint subdivision */
} geomlib_subdivision_scheme;

/* Subdivision parameters */
typedef struct geomlib_subdivision_params {
    geomlib_subdivision_scheme scheme;
    int levels;                    /* Number of subdivision levels */
    bool preserve_boundaries;      /* Keep boundary vertices fixed */
    bool project_to_limit;        /* Project to limit surface (smoother) */
} geomlib_subdivision_params;

/* Default subdivision parameters */
geomlib_subdivision_params geomlib_subdivision_params_default(
    geomlib_subdivision_scheme scheme,
    int levels
);

/* Subdivide mesh using specified scheme
 * Returns a new subdivided mesh, original is unchanged
 */
int geomlib_subdivide_mesh(
    const geomlib_trimesh *input,
    const geomlib_subdivision_params *params,
    geomlib_trimesh *output
);

/* Loop subdivision (one level)
 * Based on Charles Loop's 1987 master's thesis
 * Creates smooth surfaces from triangle meshes
 */
int geomlib_subdivide_loop(
    const geomlib_trimesh *input,
    geomlib_trimesh *output
);

/* Butterfly subdivision (one level)
 * Interpolating subdivision scheme
 * Passes through original vertices
 */
int geomlib_subdivide_butterfly(
    const geomlib_trimesh *input,
    geomlib_trimesh *output
);

/* Simple midpoint subdivision (one level)
 * Split each triangle into 4 by adding midpoints
 */
int geomlib_subdivide_midpoint(
    const geomlib_trimesh *input,
    geomlib_trimesh *output
);

/* Compute smooth limit surface position for a vertex
 * Used for projecting subdivided mesh to limit surface
 */
geomlib_vec3 geomlib_loop_limit_position(
    const geomlib_trimesh *mesh,
    uint32_t vertex_idx
);

/* Compute smooth limit normal for a vertex */
geomlib_vec3 geomlib_loop_limit_normal(
    const geomlib_trimesh *mesh,
    uint32_t vertex_idx
);

/* Helper: Compute Loop subdivision weights (beta) */
float geomlib_loop_beta(int n);  /* n = vertex valence */

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SUBDIVISION_H */
