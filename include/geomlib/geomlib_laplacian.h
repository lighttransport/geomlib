/*
 * geomlib - Cotangent Laplacian for triangle meshes
 *
 * Builds the symmetric cotangent stiffness matrix and lumped area mass matrix
 * in CSR form for PMF-based geometry processing.
 */

#ifndef GEOMLIB_LAPLACIAN_H
#define GEOMLIB_LAPLACIAN_H

#include "geomlib_sparse.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_spmat stiffness;   /* cotangent Laplacian (symmetric) */
    double *mass;              /* lumped vertex areas, vertex_count entries */
    size_t vertex_count;
} geomlib_laplacian;

geomlib_laplacian geomlib_laplacian_create(void);
void geomlib_laplacian_free(geomlib_laplacian *l);

/* Build from a triangle mesh. */
int geomlib_build_laplacian(const geomlib_trimesh *mesh, geomlib_laplacian *out);

/* Apply the operator: y = (L + scale*M) x using the lumped mass. */
int geomlib_laplacian_apply(const geomlib_laplacian *l, const double *x,
                            double scale, double *y);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_LAPLACIAN_H */
