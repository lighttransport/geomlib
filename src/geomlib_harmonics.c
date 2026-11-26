/*
 * geomlib - Manifold harmonics (see geomlib_harmonics.h)
 *
 * Clean-room implementation of the discrete Laplace-Beltrami eigenproblem
 * with cotangent stiffness and lumped area mass (Desbrun et al. 1999,
 * "Implicit fairing of irregular meshes using diffusion and curvature flow";
 * Levy, "Laplace-Beltrami eigenfunctions towards an algorithm that
 * understands geometry", 2006).
 *
 * The generalized problem L x = lambda M x is reduced to the standard
 * symmetric problem M^{-1/2} L M^{-1/2} y = lambda y and solved with the
 * existing cyclic Jacobi eigensolver. C11, wasm32 safe.
 */

#include "geomlib/geomlib_harmonics.h"
#include "geomlib/geomlib_sparse.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_harmonics_result geomlib_harmonics_result_create(size_t count,
                                                          size_t vertex_count)
{
    geomlib_harmonics_result r;
    r.eigenvalues = (double *)malloc((count ? count : 1) * sizeof(double));
    r.eigenvectors = (double *)malloc((count ? count : 1) * vertex_count * sizeof(double));
    r.count = count;
    r.vertex_count = vertex_count;
    if (!r.eigenvalues || !r.eigenvectors) geomlib_harmonics_result_free(&r);
    return r;
}

void geomlib_harmonics_result_free(geomlib_harmonics_result *r)
{
    if (!r) return;
    free(r->eigenvalues);
    free(r->eigenvectors);
    r->eigenvalues = NULL;
    r->eigenvectors = NULL;
    r->count = r->vertex_count = 0;
}

int geomlib_compute_manifold_harmonics(const geomlib_trimesh *mesh,
                                       size_t count,
                                       geomlib_harmonics_result *out)
{
    geomlib_spmat l, mass;
    size_t n, i, k;
    double *dense = NULL, *mass_diag = NULL, *eig = NULL, *vec = NULL;
    int rc;

    if (!mesh || !out || !mesh->vertices || !mesh->indices ||
        mesh->triangle_count == 0 || count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    n = mesh->vertex_count;
    if (count > n) count = n;
    if (n > 512) return GEOMLIB_MESH_ERROR_UNSUPPORTED;

    if (geomlib_cotmatrix(mesh, &l) != GEOMLIB_SPARSE_SUCCESS)
        return GEOMLIB_MESH_ERROR_INVALID;
    if (geomlib_massmatrix(mesh, &mass) != GEOMLIB_SPARSE_SUCCESS) {
        geomlib_spmat_free(&l);
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    mass_diag = (double *)calloc(n, sizeof(double));
    dense = (double *)calloc(n*n, sizeof(double));
    eig = (double *)malloc(n * sizeof(double));
    vec = (double *)malloc(n*n * sizeof(double));
    if (!mass_diag || !dense || !eig || !vec) {
        rc = GEOMLIB_MESH_ERROR_MEMORY;
        goto cleanup;
    }
    for (i = 0; i < mass.nrows; ++i) {
        for (k = mass.rowptr[i]; k < mass.rowptr[i+1]; ++k) {
            if (mass.colidx[k] == i) mass_diag[i] = mass.values[k];
        }
    }
    {
        double avg_mass = 0.0;
        for (i = 0; i < n; ++i) avg_mass += mass_diag[i];
        avg_mass /= (double)n;
        if (!(avg_mass > 0.0)) { rc = GEOMLIB_MESH_ERROR_INVALID; goto cleanup; }
        /* Isolated or degenerate vertices can have zero lumped area; give
         * them a small positive mass so the generalized problem stays
         * well-defined instead of failing the whole computation. */
        for (i = 0; i < n; ++i) {
            if (!(mass_diag[i] > 0.0)) mass_diag[i] = avg_mass * 1e-3;
        }
    }
    for (i = 0; i < l.nrows; ++i) {
        for (k = l.rowptr[i]; k < l.rowptr[i+1]; ++k) {
            size_t j = l.colidx[k];
            if (j < n) dense[i*n + j] = l.values[k] / sqrt(mass_diag[i]*mass_diag[j]);
        }
    }

    rc = geomlib_eigen_sym_dense(dense, n, eig, vec);
    if (rc != GEOMLIB_SPARSE_SUCCESS) { rc = GEOMLIB_MESH_ERROR_INVALID; goto cleanup; }

    {
        geomlib_harmonics_result result = geomlib_harmonics_result_create(count, n);
        if (!result.eigenvalues || !result.eigenvectors) {
            geomlib_harmonics_result_free(&result);
            rc = GEOMLIB_MESH_ERROR_MEMORY;
            goto cleanup;
        }
        for (i = 0; i < count; ++i) {
            result.eigenvalues[i] = eig[i];
            for (k = 0; k < n; ++k)
                result.eigenvectors[i*n + k] = vec[i*n + k] / sqrt(mass_diag[k]);
        }
        /* Release caller-owned buffers replaced by the new result. */
        free(out->eigenvalues);
        free(out->eigenvectors);
        *out = result;
    }

cleanup:
    free(mass_diag); free(dense); free(eig); free(vec);
    geomlib_spmat_free(&l);
    geomlib_spmat_free(&mass);
    return rc;
}
