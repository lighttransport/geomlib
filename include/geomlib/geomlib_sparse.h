/*
 * geomlib - Sparse linear algebra
 *
 * Compressed sparse row (CSR) matrices and iterative/direct solvers for the
 * small-to-medium sparse systems that arise in geometry processing.
 *
 * The mathematical methods are standard: conjugate gradients with Jacobi
 * preconditioning, CGLS (CG on the normal equations without forming them),
 * Householder QR least squares and the cyclic Jacobi symmetric eigensolver.
 * Clean-room implementations; no external linear algebra library.
 */

#ifndef GEOMLIB_SPARSE_H
#define GEOMLIB_SPARSE_H

#include "geomlib_mesh.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Compressed sparse row matrix (double precision). */
typedef struct {
    size_t nrows;
    size_t ncols;
    size_t nnz;
    size_t *rowptr;   /* length nrows + 1 */
    uint32_t *colidx; /* length nnz, sorted within each row */
    double *values;   /* length nnz */
} geomlib_spmat;

/* Return codes. */
#define GEOMLIB_SPARSE_SUCCESS        0
#define GEOMLIB_SPARSE_ERROR_MEMORY  (-1)
#define GEOMLIB_SPARSE_ERROR_INVALID (-2)
#define GEOMLIB_SPARSE_NO_CONVERGE   (-3)

/* Allocate an empty (zero-filled) CSR matrix. */
geomlib_spmat geomlib_spmat_create(size_t nrows, size_t ncols, size_t nnz);
void geomlib_spmat_free(geomlib_spmat *m);

/* Build a CSR matrix from raw triplets. Duplicate entries are summed and
 * column indices are sorted within each row. */
int geomlib_spmat_from_triplets(const uint32_t *rows, const uint32_t *cols,
                                const double *values, size_t count,
                                size_t nrows, size_t ncols,
                                geomlib_spmat *out);

/* y = A * x. */
int geomlib_spmat_mult_vec(const geomlib_spmat *a, const double *x, double *y);

/* out = A^T. */
int geomlib_spmat_transpose(const geomlib_spmat *a, geomlib_spmat *out);

/* Solver parameters shared by the iterative routines. */
typedef struct {
    int max_iterations;  /* default: 10 * n */
    double tolerance;    /* default: 1e-10 (relative residual) */
    bool use_jacobi;     /* default: true */
} geomlib_sparse_solve_params;

geomlib_sparse_solve_params geomlib_sparse_solve_params_default(void);

/* Solve A x = b for a square SPD matrix (Jacobi-preconditioned CG). */
int geomlib_solve_cg_sparse(const geomlib_spmat *a, const double *b, double *x,
                            const geomlib_sparse_solve_params *params);

/* Solve min ||A x - b||_2 for a possibly rectangular matrix (CGLS). */
int geomlib_solve_cgls(const geomlib_spmat *a, const double *b, double *x,
                       const geomlib_sparse_solve_params *params);

/* Dense Householder QR least squares: A is row-major m x n (m >= n). */
int geomlib_solve_qr_dense(const double *a, size_t m, size_t n,
                           const double *b, double *x);

/* Symmetric dense eigensolver (cyclic Jacobi). A is row-major n x n and is
 * not modified. Eigenvalues are ascending; eigenvectors are row-major with
 * eigenvector i stored in row i. */
int geomlib_eigen_sym_dense(const double *a, size_t n,
                            double *eigenvalues, double *eigenvectors);

/* Geometry helpers: cotangent stiffness and lumped Voronoi mass matrices
 * (double precision), with a fixed-free convention for pinned vertices. */
int geomlib_cotmatrix(const geomlib_trimesh *mesh, geomlib_spmat *out);
int geomlib_massmatrix(const geomlib_trimesh *mesh, geomlib_spmat *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SPARSE_H */
