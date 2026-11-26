/*
 * geomlib - Sparse direct solver (see geomlib_sparse_direct.h)
 *
 * Clean-room implementation: reverse Cuthill-McKee ordering followed by
 * banded Cholesky. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_sparse_direct.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static size_t build_rcm(const geomlib_spmat *a, size_t *order)
{
    size_t n = a->nrows, i, j, k = 0, head = 0, tail = 0;
    size_t *degree, *queue;
    int *visited;

    degree = (size_t *)malloc(n*sizeof(size_t));
    queue = (size_t *)malloc(n*sizeof(size_t));
    visited = (int *)calloc(n, sizeof(int));
    if (!degree || !queue || !visited) {
        free(degree); free(queue); free(visited);
        return 0;
    }
    for (i = 0; i < n; ++i) {
        degree[i] = a->rowptr[i+1] - a->rowptr[i];
        visited[i] = 0;
    }
    /* Seed at the minimum-degree vertex, then BFS with ascending-degree
     * neighbours (standard RCM variant). */
    {
        size_t seed = 0;
        for (i = 1; i < n; ++i) if (degree[i] < degree[seed]) seed = i;
        queue[tail++] = seed;
        visited[seed] = 1;
    }
    while (head < tail) {
        size_t v = queue[head++];
        order[k++] = v;
        /* Visit neighbours in natural order (no in-place column sorting: that
         * would detach colidx from values). */
        for (j = a->rowptr[v]; j < a->rowptr[v+1]; ++j) {
            size_t nb = a->colidx[j];
            if (visited[nb]) continue;
            visited[nb] = 1;
            if (tail < n) queue[tail++] = nb;
        }
    }
    for (i = 0; i < n; ++i) if (!visited[i]) order[k++] = i;
    /* Reverse for RCM. */
    for (i = 0; i < n/2; ++i) {
        size_t t = order[i];
        order[i] = order[n-1-i];
        order[n-1-i] = t;
    }
    free(degree); free(queue); free(visited);
    return n;
}

int geomlib_solve_sparse_direct(const geomlib_spmat *a, const double *b,
                                double *x, size_t *bandwidth_out)
{
    size_t n, i, j, k, bw = 0;
    size_t *order, *inv;
    double *band, *bb, *y;
    int rc = GEOMLIB_SPARSE_ERROR_INVALID;

    if (!a || !b || !x || !a->rowptr || a->nrows != a->ncols)
        return GEOMLIB_SPARSE_ERROR_INVALID;
    n = a->nrows;
    if (!n) return GEOMLIB_SPARSE_SUCCESS;
    order = (size_t *)malloc(n*sizeof(size_t));
    inv = (size_t *)malloc(n*sizeof(size_t));
    if (!order || !inv) {
        free(order); free(inv);
        return GEOMLIB_SPARSE_ERROR_MEMORY;
    }
    if (!build_rcm(a, order)) { free(order); free(inv); return GEOMLIB_SPARSE_ERROR_MEMORY; }
    for (i = 0; i < n; ++i) inv[order[i]] = i;

    /* Determine bandwidth in the reordered system. */
    for (i = 0; i < a->nrows; ++i) {
        size_t ri = inv[i];
        for (k = a->rowptr[i]; k < a->rowptr[i+1]; ++k) {
            size_t cj = inv[a->colidx[k]];
            size_t d = ri > cj ? ri-cj : cj-ri;
            if (d > bw) bw = d;
        }
    }
    if (bandwidth_out) *bandwidth_out = bw;

    band = (double *)calloc(n*(bw+1), sizeof(double));
    bb = (double *)malloc(n*sizeof(double));
    y = (double *)malloc(n*sizeof(double));
    if (!band || !bb || !y) {
        rc = GEOMLIB_SPARSE_ERROR_MEMORY;
        goto cleanup;
    }
    for (i = 0; i < n; ++i) {
        bb[i] = b[order[i]];
        for (k = a->rowptr[order[i]]; k < a->rowptr[order[i]+1]; ++k) {
            size_t cj = inv[a->colidx[k]];
            if (i >= cj && i-cj <= bw)
                band[i*(bw+1) + (i-cj)] = a->values[k];
        }
    }

    /* Banded Cholesky: L[i*(bw+1)+j] holds L[i][i-j]. */
    for (i = 0; i < n; ++i) {
        size_t lo = i > bw ? i-bw : 0;
        for (j = lo; j < i; ++j) {
            double sum = band[i*(bw+1)+(i-j)];
            size_t q;
            /* L[i][j] = (A[i][j] - sum_{k<j} L[i][k] L[j][k]) / L[j][j].
             * Both factors live in the stored lower band: k >= lo = i-bw
             * guarantees i-k <= bw and (j-k) < bw since j <= i-1. */
            for (q = lo; q < j; ++q)
                sum -= band[i*(bw+1)+(i-q)] * band[j*(bw+1)+(j-q)];
            if (fabs(band[j*(bw+1)]) < 1e-14) { rc = GEOMLIB_SPARSE_ERROR_INVALID; goto cleanup; }
            band[i*(bw+1)+(i-j)] = sum / band[j*(bw+1)];
        }
        {
            double d = band[i*(bw+1)];
            size_t q;
            for (q = lo; q < i; ++q) {
                size_t iq = i-q;
                if (iq <= bw) d -= band[i*(bw+1)+iq]*band[i*(bw+1)+iq];
            }
            if (!(d > 1e-14)) { rc = GEOMLIB_SPARSE_ERROR_INVALID; goto cleanup; }
            band[i*(bw+1)] = sqrt(d);
        }
    }

    for (i = 0; i < n; ++i) {
        double s = bb[i];
        size_t lo = i > bw ? i-bw : 0;
        for (j = lo; j < i; ++j) {
            size_t ij = i-j;
            if (ij <= bw) s -= band[i*(bw+1)+ij]*y[j];
        }
        y[i] = s / band[i*(bw+1)];
    }
    for (i = n; i-- > 0;) {
        double s = y[i];
        size_t hi = i+bw < n ? i+bw : n-1;
        for (j = i+1; j <= hi; ++j) {
            size_t ji = j-i;
            if (ji <= bw) s -= band[j*(bw+1)+ji]*x[order[j]];
        }
        x[order[i]] = s / band[i*(bw+1)];
    }
    rc = GEOMLIB_SPARSE_SUCCESS;

cleanup:
    free(band); free(bb); free(y);
    free(order); free(inv);
    return rc;
}
