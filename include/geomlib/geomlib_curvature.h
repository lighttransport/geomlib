/*
 * geomlib - Curvature Computation
 * Discrete curvature estimation on triangle meshes
 * Based on Meyer et al. "Discrete Differential-Geometry Operators for Triangulated 2-Manifolds" (2003)
 */

#ifndef GEOMLIB_CURVATURE_H
#define GEOMLIB_CURVATURE_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Curvature values at a vertex */
typedef struct geomlib_vertex_curvature {
    float mean;              /* Mean curvature H */
    float gaussian;          /* Gaussian curvature K */
    float principal_max;     /* Maximum principal curvature k1 */
    float principal_min;     /* Minimum principal curvature k2 */
    geomlib_vec3 principal_dir_max;  /* Direction of k1 */
    geomlib_vec3 principal_dir_min;  /* Direction of k2 */
} geomlib_vertex_curvature;

/* Curvature data for entire mesh */
typedef struct geomlib_mesh_curvature {
    geomlib_vertex_curvature *curvatures;
    size_t vertex_count;
    float *vertex_areas;     /* Voronoi areas for each vertex */
} geomlib_mesh_curvature;

/* Create curvature structure */
geomlib_mesh_curvature geomlib_curvature_create(size_t vertex_count);

/* Free curvature structure */
void geomlib_curvature_free(geomlib_mesh_curvature *curv);

/* Compute all curvatures for mesh */
int geomlib_compute_curvature(
    const geomlib_trimesh *mesh,
    geomlib_mesh_curvature *curvature
);

/* Compute mean curvature only (faster) */
int geomlib_compute_mean_curvature(
    const geomlib_trimesh *mesh,
    float *mean_curvatures
);

/* Compute Gaussian curvature only */
int geomlib_compute_gaussian_curvature(
    const geomlib_trimesh *mesh,
    float *gaussian_curvatures
);

/* Compute principal curvatures from mean and Gaussian */
void geomlib_compute_principal_curvatures(
    float mean,
    float gaussian,
    float *k1,
    float *k2
);

/* Helper: Compute mixed Voronoi area for vertex */
float geomlib_compute_voronoi_area(
    const geomlib_trimesh *mesh,
    uint32_t vertex_idx
);

/* Helper: Compute cotan weight for edge */
float geomlib_cotan_weight(
    const geomlib_vec3 *v0,
    const geomlib_vec3 *v1,
    const geomlib_vec3 *v2
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CURVATURE_H */
