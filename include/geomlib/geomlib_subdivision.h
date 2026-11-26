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
    GEOMLIB_SUBDIV_CATMULL_CLARK, /* Catmull-Clark (quadrilateral refinement) */
    GEOMLIB_SUBDIV_MIDPOINT       /* Simple midpoint subdivision */
} geomlib_subdivision_scheme;

/* Subdivision parameters */
typedef struct geomlib_subdivision_params {
    geomlib_subdivision_scheme scheme;
    int levels;                    /* Number of subdivision levels */
    bool preserve_boundaries;      /* Loop: freeze boundary vertices */
    bool project_to_limit;         /* Loop: project interior vertices to the limit surface */
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

/* Modified Butterfly subdivision (Zorin et al.), one level. Interpolating:
 * control vertices are unchanged. Regular edges use the 8-point stencil
 * 1/2(a+b) + 1/8(c+d) - 1/16(sum of four outer vertices); edges touching an
 * extraordinary vertex use Zorin's extraordinary masks, and boundary edges use
 * the linear midpoint. */
int geomlib_subdivide_butterfly(
    const geomlib_trimesh *input,
    geomlib_trimesh *output
);

/* Catmull-Clark subdivision (one level). Faces are triangulated after the
 * quad split; boundary vertices use the standard crease rule. */
int geomlib_subdivide_catmull_clark(
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

/* Loop limit position using the evaluation mask
 * alpha_n = (3/8 + cos(2*pi/n)/4)^2 + 3/8 and epsilon_n = 3/(11 - 8 alpha_n):
 * v_inf = epsilon_n v + (1 - epsilon_n) * mean(one-ring). Boundary vertices keep
 * their position. */
geomlib_vec3 geomlib_loop_limit_position(
    const geomlib_trimesh *mesh,
    uint32_t vertex_idx
);

/* Loop limit normal from the tangent masks tau_k = (cos(2*pi*k/n), sin(2*pi*k/n))
 * over the cyclically ordered one-ring. Boundary vertices use the area-weighted
 * vertex normal. */
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
