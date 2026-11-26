/*
 * geomlib - Manifold harmonics (Laplace-Beltrami eigenfunctions)
 *
 * Computes the lowest eigenvalues and eigenfunctions of the cotangent
 * Laplace-Beltrami operator with lumped mass matrix, solved as a generalized
 * symmetric eigenproblem L x = lambda M x.
 */

#ifndef GEOMLIB_HARMONICS_H
#define GEOMLIB_HARMONICS_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double *eigenvalues;   /* Ascending, length count */
    double *eigenvectors;  /* Row-major: eigenvector i in row i, vertex count columns */
    size_t count;
    size_t vertex_count;
} geomlib_harmonics_result;

geomlib_harmonics_result geomlib_harmonics_result_create(size_t count,
                                                          size_t vertex_count);
void geomlib_harmonics_result_free(geomlib_harmonics_result *r);

/*
 * Compute `count` lowest Laplace-Beltrami eigenpairs of a triangle mesh.
 * Uses the cotangent stiffness matrix and lumped area mass matrix.
 */
int geomlib_compute_manifold_harmonics(const geomlib_trimesh *mesh,
                                       size_t count,
                                       geomlib_harmonics_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_HARMONICS_H */
