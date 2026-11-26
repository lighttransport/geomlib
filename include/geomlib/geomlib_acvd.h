/*
 * geomlib - ACVD-style remeshing
 *
 * Approximated Centroidal Voronoi Diagram remeshing: relax a set of surface
 * sites toward the centroids of their restricted Voronoi cells (Lloyd
 * iterations) and triangulate the relaxed sites in local tangent planes.
 */

#ifndef GEOMLIB_ACVD_H
#define GEOMLIB_ACVD_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t target_vertices;  /* Target number of vertices (sites) */
    int lloyd_iterations;      /* Default 5 */
    uint32_t seed;             /* RNG seed (0 = deterministic default) */
} geomlib_acvd_params;

geomlib_acvd_params geomlib_acvd_params_default(uint32_t target_vertices);

/*
 * Remesh a closed surface with approximately `target_vertices` vertices.
 * The output topology is rebuilt; geometry is kept on the input surface.
 */
int geomlib_remesh_acvd(const geomlib_trimesh *input,
                        const geomlib_acvd_params *params,
                        geomlib_trimesh *output);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ACVD_H */
