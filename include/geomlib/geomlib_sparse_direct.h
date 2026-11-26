/*
 * geomlib - Sparse direct solver for symmetric positive-definite systems
 *
 * Uses reverse Cuthill-McKee bandwidth reduction followed by banded
 * Cholesky factorization and triangular solves.
 */

#ifndef GEOMLIB_SPARSE_DIRECT_H
#define GEOMLIB_SPARSE_DIRECT_H

#include "geomlib_sparse.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Solve A x = b for a square SPD CSR matrix. Returns GEOMLIB_SPARSE_SUCCESS
 * on success. `bandwidth_out` optionally receives the RCM bandwidth. */
int geomlib_solve_sparse_direct(const geomlib_spmat *a, const double *b,
                                double *x, size_t *bandwidth_out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SPARSE_DIRECT_H */
