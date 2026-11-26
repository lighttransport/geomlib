/*
 * geomlib - Implicit mean-curvature flow (see geomlib_fairing.h)
 *
 * Clean-room implementation of Desbrun et al., "Implicit Fairing of Irregular
 * Meshes using Diffusion and Curvature Flow", SIGGRAPH 1999. Each step solves
 * (M + lambda L) x = M x_old with the sparse CG solver.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_fairing.h"
#include "geomlib/geomlib_sparse.h"
#include "geomlib/geomlib_topology.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_fairing_params geomlib_fairing_params_default(void)
{
    geomlib_fairing_params p;
    p.lambda = 0.1f;
    p.iterations = 5;
    p.preserve_boundaries = true;
    return p;
}

/* Build A = M + lambda*L, with fixed rows for boundary vertices. */
static int build_matrix(const geomlib_trimesh *mesh, float lambda,
                        const bool *fixed, const double *mass_diag,
                        geomlib_spmat *out)
{
    geomlib_spmat l;
    size_t *rows, *cols, cap, count = 0;
    double *vals;
    size_t r, k;

    if (geomlib_cotmatrix(mesh, &l) != GEOMLIB_SPARSE_SUCCESS)
        return GEOMLIB_MESH_ERROR_INVALID;
    cap = l.nnz + mesh->vertex_count*2 + 1;
    rows = (size_t *)malloc(cap*sizeof(size_t));
    cols = (size_t *)malloc(cap*sizeof(size_t));
    vals = (double *)malloc(cap*sizeof(double));
    if (!rows || !cols || !vals) {
        free(rows); free(cols); free(vals);
        geomlib_spmat_free(&l);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (r = 0; r < mesh->vertex_count; ++r) {
        if (fixed && fixed[r]) {
            rows[count]=r; cols[count]=r; vals[count++]=1.0;
            continue;
        }
        rows[count]=r; cols[count]=r; vals[count++]=mass_diag[r];
        for (k = l.rowptr[r]; k < l.rowptr[r+1]; ++k) {
            size_t c = l.colidx[k];
            if (fixed && fixed[c]) continue;
            rows[count]=r; cols[count]=c; vals[count++]=lambda*l.values[k];
        }
    }
    {
        uint32_t *rr = (uint32_t *)malloc(count*sizeof(uint32_t));
        uint32_t *cc = (uint32_t *)malloc(count*sizeof(uint32_t));
        if (!rr || !cc) {
            free(rr); free(cc); free(rows); free(cols); free(vals);
            geomlib_spmat_free(&l);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        for (k = 0; k < count; ++k) { rr[k]=(uint32_t)rows[k]; cc[k]=(uint32_t)cols[k]; }
        {
            int rc = geomlib_spmat_from_triplets(rr, cc, vals, count,
                                                 mesh->vertex_count,
                                                 mesh->vertex_count, out);
            (void)rc;
        }
        free(rr); free(cc);
    }
    free(rows); free(cols); free(vals);
    geomlib_spmat_free(&l);
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_fair_mesh(geomlib_trimesh *mesh,
                      const geomlib_fairing_params *params)
{
    size_t V, i;
    int iter, iters;
    float lambda;
    bool *fixed = NULL;
    double *mass_diag = NULL, *rhs = NULL, *sol = NULL;
    geomlib_spmat mass;
    geomlib_sparse_solve_params sp;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!mesh || !mesh->vertices || !mesh->indices || mesh->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    V = mesh->vertex_count;
    lambda = params ? params->lambda : 0.1f;
    if (!(lambda > 0.0f)) lambda = 0.1f;
    iters = params ? params->iterations : 5;
    if (iters < 0) iters = 0;

    if (params && params->preserve_boundaries) {
        geomlib_boundary_loops loops = geomlib_extract_boundary_loops(mesh);
        fixed = (bool *)calloc(V, sizeof(bool));
        if (!fixed) { geomlib_boundary_loops_free(&loops); return GEOMLIB_MESH_ERROR_MEMORY; }
        for (i = 0; i < loops.loop_count; ++i) {
            const geomlib_boundary_loop *loop = &loops.loops[i];
            size_t j;
            for (j = 0; j < loop->vertex_count; ++j)
                if (loop->vertices[j] < V) fixed[loop->vertices[j]] = true;
        }
        geomlib_boundary_loops_free(&loops);
    }

    if (geomlib_massmatrix(mesh, &mass) != GEOMLIB_SPARSE_SUCCESS) {
        free(fixed);
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    mass_diag = (double *)calloc(V, sizeof(double));
    rhs = (double *)calloc(V, sizeof(double));
    sol = (double *)calloc(V, sizeof(double));
    if (!mass_diag || !rhs || !sol) {
        rc = GEOMLIB_MESH_ERROR_MEMORY;
        goto cleanup;
    }
    for (i = 0; i < mass.nrows; ++i) {
        for (size_t k = mass.rowptr[i]; k < mass.rowptr[i+1]; ++k)
            if (mass.colidx[k] == i) mass_diag[i] = mass.values[k];
        if (!(mass_diag[i] > 0.0)) mass_diag[i] = DBL_EPSILON;
    }

    sp = geomlib_sparse_solve_params_default();
    sp.max_iterations = 4000;
    sp.tolerance = 1e-6f;

    for (iter = 0; iter < iters; ++iter) {
        geomlib_spmat A = {0};
        for (i = 0; i < V; ++i) {
            rhs[i] = mass_diag[i] * mesh->vertices[i].x;
            if (fixed && fixed[i]) rhs[i] = mesh->vertices[i].x;
        }
        rc = build_matrix(mesh, lambda, fixed, mass_diag, &A);
        if (rc != GEOMLIB_MESH_SUCCESS) goto cleanup;
        if (geomlib_solve_cg_sparse(&A, rhs, sol, &sp) != GEOMLIB_SPARSE_SUCCESS) {
            geomlib_spmat_free(&A);
            rc = GEOMLIB_MESH_ERROR_INVALID;
            goto cleanup;
        }
        geomlib_spmat_free(&A);
        for (i = 0; i < V; ++i) mesh->vertices[i].x = (float)sol[i];

        for (i = 0; i < V; ++i) {
            rhs[i] = mass_diag[i] * mesh->vertices[i].y;
            if (fixed && fixed[i]) rhs[i] = mesh->vertices[i].y;
        }
        rc = build_matrix(mesh, lambda, fixed, mass_diag, &A);
        if (rc != GEOMLIB_MESH_SUCCESS) goto cleanup;
        if (geomlib_solve_cg_sparse(&A, rhs, sol, &sp) != GEOMLIB_SPARSE_SUCCESS) {
            geomlib_spmat_free(&A);
            rc = GEOMLIB_MESH_ERROR_INVALID;
            goto cleanup;
        }
        geomlib_spmat_free(&A);
        for (i = 0; i < V; ++i) mesh->vertices[i].y = (float)sol[i];

        for (i = 0; i < V; ++i) {
            rhs[i] = mass_diag[i] * mesh->vertices[i].z;
            if (fixed && fixed[i]) rhs[i] = mesh->vertices[i].z;
        }
        rc = build_matrix(mesh, lambda, fixed, mass_diag, &A);
        if (rc != GEOMLIB_MESH_SUCCESS) goto cleanup;
        if (geomlib_solve_cg_sparse(&A, rhs, sol, &sp) != GEOMLIB_SPARSE_SUCCESS) {
            geomlib_spmat_free(&A);
            rc = GEOMLIB_MESH_ERROR_INVALID;
            goto cleanup;
        }
        geomlib_spmat_free(&A);
        for (i = 0; i < V; ++i) mesh->vertices[i].z = (float)sol[i];
    }
    rc = GEOMLIB_MESH_SUCCESS;

cleanup:
    free(fixed); free(mass_diag); free(rhs); free(sol);
    geomlib_spmat_free(&mass);
    return rc;
}
