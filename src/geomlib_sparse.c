/*
 * geomlib - Sparse linear algebra (see geomlib_sparse.h)
 *
 * Clean-room implementations of standard numerical methods. No external
 * dependencies, no long double, C11 + libm only (wasm32 compatible).
 */

#include "geomlib/geomlib_sparse.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------- */
/* CSR matrix                                                                */
/* ------------------------------------------------------------------------- */

geomlib_spmat geomlib_spmat_create(size_t nrows, size_t ncols, size_t nnz)
{
    geomlib_spmat m;
    m.nrows = nrows;
    m.ncols = ncols;
    m.nnz = nnz;
    m.rowptr = (size_t *)calloc(nrows + 1, sizeof(size_t));
    m.colidx = (uint32_t *)calloc(nnz ? nnz : 1, sizeof(uint32_t));
    m.values = (double *)calloc(nnz ? nnz : 1, sizeof(double));
    if (!m.rowptr || !m.colidx || !m.values) {
        free(m.rowptr); free(m.colidx); free(m.values);
        m.rowptr = NULL; m.colidx = NULL; m.values = NULL;
        m.nrows = m.ncols = m.nnz = 0;
    }
    return m;
}

void geomlib_spmat_free(geomlib_spmat *m)
{
    if (!m) return;
    free(m->rowptr);
    free(m->colidx);
    free(m->values);
    m->rowptr = NULL;
    m->colidx = NULL;
    m->values = NULL;
    m->nrows = m->ncols = m->nnz = 0;
}

typedef struct {
    uint32_t row;
    uint32_t col;
    double val;
} triplet_entry;

static int compare_triplet(const void *pa, const void *pb)
{
    const triplet_entry *a = (const triplet_entry *)pa;
    const triplet_entry *b = (const triplet_entry *)pb;
    if (a->row < b->row) return -1;
    if (a->row > b->row) return 1;
    if (a->col < b->col) return -1;
    if (a->col > b->col) return 1;
    return 0;
}

int geomlib_spmat_from_triplets(const uint32_t *rows, const uint32_t *cols,
                                const double *values, size_t count,
                                size_t nrows, size_t ncols,
                                geomlib_spmat *out)
{
    triplet_entry *entries;
    size_t i, write = 0, next_row = 0;
    uint32_t lastrow = 0, lastcol = 0;
    bool have_last = false;

    if (!out || !rows || !cols || !values) return GEOMLIB_SPARSE_ERROR_INVALID;
    if (nrows == 0 || (count > 0 && ncols == 0)) return GEOMLIB_SPARSE_ERROR_INVALID;
    out->nrows = nrows;
    out->ncols = ncols;
    out->rowptr = NULL;
    out->colidx = NULL;
    out->values = NULL;
    out->nnz = 0;

    for (i = 0; i < count; ++i) {
        if (rows[i] >= nrows || (ncols > 0 && cols[i] >= ncols)) {
            return GEOMLIB_SPARSE_ERROR_INVALID;
        }
    }

    entries = (triplet_entry *)malloc((count ? count : 1) * sizeof(triplet_entry));
    if (!entries) return GEOMLIB_SPARSE_ERROR_MEMORY;
    for (i = 0; i < count; ++i) {
        entries[i].row = rows[i];
        entries[i].col = cols[i];
        entries[i].val = values[i];
    }
    qsort(entries, count, sizeof(triplet_entry), compare_triplet);

    out->rowptr = (size_t *)calloc(nrows + 1, sizeof(size_t));
    out->colidx = (uint32_t *)malloc((count ? count : 1) * sizeof(uint32_t));
    out->values = (double *)malloc((count ? count : 1) * sizeof(double));
    if (!out->rowptr || !out->colidx || !out->values) {
        free(entries);
        geomlib_spmat_free(out);
        return GEOMLIB_SPARSE_ERROR_MEMORY;
    }

    for (i = 0; i < count; ++i) {
        uint32_t r = entries[i].row;
        uint32_t c = entries[i].col;
        if (!have_last || r != lastrow) {
            while (next_row <= r) {
                out->rowptr[next_row] = write;
                ++next_row;
            }
        }
        if (have_last && lastrow == r && lastcol == c) {
            out->values[write - 1] += entries[i].val;
        } else {
            out->colidx[write] = c;
            out->values[write] = entries[i].val;
            ++write;
            lastrow = r;
            lastcol = c;
            have_last = true;
        }
    }
    while (next_row <= nrows) {
        out->rowptr[next_row] = write;
        ++next_row;
    }
    out->nnz = write;

    free(entries);
    return GEOMLIB_SPARSE_SUCCESS;
}

int geomlib_spmat_mult_vec(const geomlib_spmat *a, const double *x, double *y)
{
    size_t r, k;
    if (!a || !x || !y || !a->rowptr) return GEOMLIB_SPARSE_ERROR_INVALID;
    for (r = 0; r < a->nrows; ++r) {
        double sum = 0.0;
        for (k = a->rowptr[r]; k < a->rowptr[r + 1]; ++k) {
            sum += a->values[k] * x[a->colidx[k]];
        }
        y[r] = sum;
    }
    return GEOMLIB_SPARSE_SUCCESS;
}

int geomlib_spmat_transpose(const geomlib_spmat *a, geomlib_spmat *out)
{
    size_t *counts, *next;
    size_t i, k;
    if (!a || !out || !a->rowptr) return GEOMLIB_SPARSE_ERROR_INVALID;

    counts = (size_t *)calloc(a->ncols + 1, sizeof(size_t));
    if (!counts) return GEOMLIB_SPARSE_ERROR_MEMORY;
    for (i = 0; i < a->nnz; ++i) counts[a->colidx[i] + 1]++;
    for (i = 1; i <= a->ncols; ++i) counts[i] += counts[i - 1];

    out->nrows = a->ncols;
    out->ncols = a->nrows;
    out->nnz = a->nnz;
    out->rowptr = counts;
    out->colidx = (uint32_t *)malloc((a->nnz ? a->nnz : 1) * sizeof(uint32_t));
    out->values = (double *)malloc((a->nnz ? a->nnz : 1) * sizeof(double));
    if (!out->colidx || !out->values) {
        geomlib_spmat_free(out);
        return GEOMLIB_SPARSE_ERROR_MEMORY;
    }

    next = (size_t *)malloc((a->ncols + 1) * sizeof(size_t));
    if (!next) {
        geomlib_spmat_free(out);
        return GEOMLIB_SPARSE_ERROR_MEMORY;
    }
    memcpy(next, counts, (a->ncols + 1) * sizeof(size_t));
    for (i = 0; i < a->nrows; ++i) {
        for (k = a->rowptr[i]; k < a->rowptr[i + 1]; ++k) {
            size_t dest = next[a->colidx[k]]++;
            out->colidx[dest] = (uint32_t)i;
            out->values[dest] = a->values[k];
        }
    }
    free(next);
    return GEOMLIB_SPARSE_SUCCESS;
}

/* ------------------------------------------------------------------------- */
/* Iterative solvers                                                         */
/* ------------------------------------------------------------------------- */

geomlib_sparse_solve_params geomlib_sparse_solve_params_default(void)
{
    geomlib_sparse_solve_params p;
    p.max_iterations = 0;
    p.tolerance = 1e-10;
    p.use_jacobi = true;
    return p;
}

static double vec_dot(const double *a, const double *b, size_t n)
{
    double s = 0.0;
    size_t i;
    for (i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

static double vec_norm(const double *a, size_t n)
{
    return sqrt(vec_dot(a, a, n));
}

static void apply_jacobi(const geomlib_spmat *a, const double *r, double *z,
                         bool enabled)
{
    size_t i, k;
    if (!enabled) {
        memcpy(z, r, a->nrows * sizeof(double));
        return;
    }
    for (i = 0; i < a->nrows; ++i) {
        double d = 0.0;
        for (k = a->rowptr[i]; k < a->rowptr[i + 1]; ++k) {
            if (a->colidx[k] == i) d = a->values[k];
        }
        z[i] = (fabs(d) > DBL_EPSILON) ? r[i] / d : r[i];
    }
}

int geomlib_solve_cg_sparse(const geomlib_spmat *a, const double *b, double *x,
                            const geomlib_sparse_solve_params *params)
{
    size_t n, i, iter;
    int maxit;
    double tol, bnorm, rz, rzold, alpha, beta, residual;
    bool jacobi;
    double *r, *z, *p, *ap;

    if (!a || !b || !x || !a->rowptr || a->nrows != a->ncols) {
        return GEOMLIB_SPARSE_ERROR_INVALID;
    }
    n = a->nrows;
    maxit = params && params->max_iterations > 0 ? params->max_iterations : (int)(10 * n + 10);
    tol = params ? params->tolerance : 1e-10;
    jacobi = params ? params->use_jacobi : true;

    /* calloc keeps cross-TU static analysis happy: the buffers are defined
     * before any external call is allowed to write into them. */
    r = (double *)calloc(n ? n : 1, sizeof(double));
    z = (double *)calloc(n ? n : 1, sizeof(double));
    p = (double *)calloc(n ? n : 1, sizeof(double));
    ap = (double *)calloc(n ? n : 1, sizeof(double));
    if (!r || !z || !p || !ap) {
        free(r); free(z); free(p); free(ap);
        return GEOMLIB_SPARSE_ERROR_MEMORY;
    }

    geomlib_spmat_mult_vec(a, x, ap);
    for (i = 0; i < n; ++i) r[i] = b[i] - ap[i];
    bnorm = vec_norm(b, n);
    if (bnorm == 0.0) bnorm = 1.0;

    apply_jacobi(a, r, z, jacobi);
    memcpy(p, z, n * sizeof(double));
    rz = vec_dot(r, z, n);

    for (iter = 0; iter < (size_t)maxit; ++iter) {
        double denom;
        if ((residual = vec_norm(r, n)) <= tol * bnorm) break;
        geomlib_spmat_mult_vec(a, p, ap);
        denom = vec_dot(p, ap, n);
        if (fabs(denom) < DBL_EPSILON) break;
        alpha = rz / denom;
        for (i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        apply_jacobi(a, r, z, jacobi);
        rzold = rz;
        rz = vec_dot(r, z, n);
        if (fabs(rzold) < DBL_EPSILON) break;
        beta = rz / rzold;
        for (i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
    }
    residual = vec_norm(r, n);
    free(r); free(z); free(p); free(ap);
    return residual <= tol * bnorm ? GEOMLIB_SPARSE_SUCCESS : GEOMLIB_SPARSE_NO_CONVERGE;
}

int geomlib_solve_cgls(const geomlib_spmat *a, const double *b, double *x,
                       const geomlib_sparse_solve_params *params)
{
    geomlib_spmat at;
    size_t m, n, i, iter;
    int maxit;
    double tol, snorm, snorm0, alpha, beta, qnorm;
    double *r, *s, *p, *q;
    int rc;

    if (!a || !b || !x || !a->rowptr) return GEOMLIB_SPARSE_ERROR_INVALID;
    m = a->nrows;
    n = a->ncols;
    rc = geomlib_spmat_transpose(a, &at);
    if (rc != GEOMLIB_SPARSE_SUCCESS) return rc;

    maxit = params && params->max_iterations > 0 ? params->max_iterations : (int)(10 * n + 10);
    tol = params ? params->tolerance : 1e-10;

    r = (double *)calloc(m ? m : 1, sizeof(double));
    s = (double *)calloc(n ? n : 1, sizeof(double));
    p = (double *)calloc(n ? n : 1, sizeof(double));
    q = (double *)calloc(m ? m : 1, sizeof(double));
    if (!r || !s || !p || !q) {
        free(r); free(s); free(p); free(q);
        geomlib_spmat_free(&at);
        return GEOMLIB_SPARSE_ERROR_MEMORY;
    }

    geomlib_spmat_mult_vec(a, x, q);
    for (i = 0; i < m; ++i) r[i] = b[i] - q[i];

    geomlib_spmat_mult_vec(&at, r, s);
    snorm0 = vec_norm(s, n);
    if (snorm0 == 0.0) snorm0 = 1.0;
    memcpy(p, s, n * sizeof(double));
    snorm = vec_dot(s, s, n);

    for (iter = 0; iter < (size_t)maxit; ++iter) {
        double snorm_old;
        if (sqrt(snorm) <= tol * snorm0) break;
        geomlib_spmat_mult_vec(a, p, q);
        qnorm = vec_dot(q, q, m);
        if (qnorm < DBL_EPSILON) break;
        alpha = snorm / qnorm;
        for (i = 0; i < n; ++i) x[i] += alpha * p[i];
        for (i = 0; i < m; ++i) r[i] -= alpha * q[i];
        geomlib_spmat_mult_vec(&at, r, s);
        snorm_old = snorm;
        snorm = vec_dot(s, s, n);
        if (fabs(snorm_old) < DBL_EPSILON) break;
        beta = snorm / snorm_old;
        for (i = 0; i < n; ++i) p[i] = s[i] + beta * p[i];
    }
    free(r); free(s); free(p); free(q);
    geomlib_spmat_free(&at);
    return sqrt(snorm) <= tol * snorm0 ? GEOMLIB_SPARSE_SUCCESS
                                       : GEOMLIB_SPARSE_NO_CONVERGE;
}

/* ------------------------------------------------------------------------- */
/* Dense direct solvers                                                      */
/* ------------------------------------------------------------------------- */

int geomlib_solve_qr_dense(const double *a, size_t m, size_t n,
                           const double *b, double *x)
{
    double *r, *y;
    size_t i, j, k;

    if (!a || !b || !x || m < n || n == 0) return GEOMLIB_SPARSE_ERROR_INVALID;
    r = (double *)malloc(m * n * sizeof(double));
    y = (double *)malloc(m * sizeof(double));
    if (!r || !y) {
        free(r); free(y);
        return GEOMLIB_SPARSE_ERROR_MEMORY;
    }
    memcpy(r, a, m * n * sizeof(double));
    memcpy(y, b, m * sizeof(double));

    for (k = 0; k < n; ++k) {
        double normx = 0.0, alpha, vnorm2 = 0.0, factor;
        double *v;
        for (i = k; i < m; ++i) normx += r[i * n + k] * r[i * n + k];
        normx = sqrt(normx);
        if (normx < DBL_EPSILON) continue;
        alpha = (r[k * n + k] >= 0.0) ? -normx : normx;
        v = (double *)malloc((m - k) * sizeof(double));
        if (!v) { free(r); free(y); return GEOMLIB_SPARSE_ERROR_MEMORY; }
        for (i = k; i < m; ++i) {
            v[i - k] = r[i * n + k];
            if (i == k) v[i - k] -= alpha;
            vnorm2 += v[i - k] * v[i - k];
        }
        if (vnorm2 > DBL_EPSILON) {
            for (j = k; j < n; ++j) {
                double dot = 0.0;
                for (i = k; i < m; ++i) dot += v[i - k] * r[i * n + j];
                factor = 2.0 * dot / vnorm2;
                for (i = k; i < m; ++i) r[i * n + j] -= factor * v[i - k];
            }
            {
                double dot = 0.0;
                for (i = k; i < m; ++i) dot += v[i - k] * y[i];
                factor = 2.0 * dot / vnorm2;
                for (i = k; i < m; ++i) y[i] -= factor * v[i - k];
            }
        }
        free(v);
    }

    for (i = n; i-- > 0;) {
        double sum = y[i];
        for (j = i + 1; j < n; ++j) sum -= r[i * n + j] * x[j];
        if (fabs(r[i * n + i]) < DBL_EPSILON) {
            free(r); free(y);
            return GEOMLIB_SPARSE_ERROR_INVALID;
        }
        x[i] = sum / r[i * n + i];
    }

    free(r); free(y);
    return GEOMLIB_SPARSE_SUCCESS;
}

int geomlib_eigen_sym_dense(const double *a, size_t n,
                            double *eigenvalues, double *eigenvectors)
{
    double *work, *v;
    size_t i, j, p = 0, q = 0, sweep;
    const size_t max_sweeps = 100 * n + 100;

    if (!a || !eigenvalues || !eigenvectors || n == 0) {
        return GEOMLIB_SPARSE_ERROR_INVALID;
    }
    work = (double *)malloc(n * n * sizeof(double));
    v = (double *)malloc(n * n * sizeof(double));
    if (!work || !v) {
        free(work); free(v);
        return GEOMLIB_SPARSE_ERROR_MEMORY;
    }
    memcpy(work, a, n * n * sizeof(double));
    memset(v, 0, n * n * sizeof(double));
    for (i = 0; i < n; ++i) v[i * n + i] = 1.0;
    memset(eigenvectors, 0, n * n * sizeof(double));

    for (sweep = 0; sweep < max_sweeps; ++sweep) {
        double maxoff = 0.0;
        double theta, t, c, s, app, aqq, apq;
        for (i = 0; i < n; ++i) {
            for (j = i + 1; j < n; ++j) {
                double val = fabs(work[i * n + j]);
                if (val > maxoff) { maxoff = val; p = i; q = j; }
            }
        }
        if (maxoff < 1e-15) break;
        app = work[p * n + p];
        aqq = work[q * n + q];
        apq = work[p * n + q];
        theta = (aqq - app) / (2.0 * apq);
        if (theta >= 0.0) t = 1.0 / (theta + sqrt(1.0 + theta * theta));
        else t = -1.0 / (-theta + sqrt(1.0 + theta * theta));
        c = 1.0 / sqrt(1.0 + t * t);
        s = t * c;

        work[p * n + p] = app - t * apq;
        work[q * n + q] = aqq + t * apq;
        work[p * n + q] = work[q * n + p] = 0.0;
        for (i = 0; i < n; ++i) {
            if (i == p || i == q) continue;
            double aip = work[i * n + p];
            double aiq = work[i * n + q];
            work[i * n + p] = work[p * n + i] = c * aip - s * aiq;
            work[i * n + q] = work[q * n + i] = s * aip + c * aiq;
        }
        for (i = 0; i < n; ++i) {
            double vip = v[i * n + p];
            double viq = v[i * n + q];
            v[i * n + p] = c * vip - s * viq;
            v[i * n + q] = s * vip + c * viq;
        }
    }

    for (i = 0; i < n; ++i) eigenvalues[i] = work[i * n + i];
    /* Sort ascending; swap eigenvector columns alongside eigenvalues. */
    for (i = 0; i < n; ++i) {
        size_t best = i;
        for (j = i + 1; j < n; ++j) {
            if (eigenvalues[j] < eigenvalues[best]) best = j;
        }
        if (best != i) {
            double tmp = eigenvalues[i];
            eigenvalues[i] = eigenvalues[best];
            eigenvalues[best] = tmp;
            for (j = 0; j < n; ++j) {
                tmp = v[j * n + i];
                v[j * n + i] = v[j * n + best];
                v[j * n + best] = tmp;
            }
        }
    }
    /* Return eigenvectors as rows: row i is eigenvector i. */
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) eigenvectors[i * n + j] = v[j * n + i];
    }

    free(work); free(v);
    return GEOMLIB_SPARSE_SUCCESS;
}

/* ------------------------------------------------------------------------- */
/* Geometry operators                                                        */
/* ------------------------------------------------------------------------- */

static double cotan_at(const geomlib_vec3 *a, const geomlib_vec3 *b,
                       const geomlib_vec3 *c)
{
    /* Cotangent of the angle at vertex b in triangle (a, b, c). */
    geomlib_vec3 ba = {a->x - b->x, a->y - b->y, a->z - b->z};
    geomlib_vec3 bc = {c->x - b->x, c->y - b->y, c->z - b->z};
    double dot = (double)ba.x * bc.x + (double)ba.y * bc.y + (double)ba.z * bc.z;
    geomlib_vec3 cr = {
        (double)ba.y * bc.z - (double)ba.z * bc.y,
        (double)ba.z * bc.x - (double)ba.x * bc.z,
        (double)ba.x * bc.y - (double)ba.y * bc.x
    };
    double cross = sqrt((double)cr.x * cr.x + (double)cr.y * cr.y + (double)cr.z * cr.z);
    if (fabs(cross) < 1e-30) return 0.0;
    return dot / cross;
}

int geomlib_cotmatrix(const geomlib_trimesh *mesh, geomlib_spmat *out)
{
    size_t tri, ntrip = 0, cap, n;
    uint32_t *rows, *cols;
    double *vals;
    int rc;

    if (!mesh || !out || !mesh->indices) return GEOMLIB_SPARSE_ERROR_INVALID;
    n = mesh->vertex_count;
    cap = mesh->triangle_count * 9 + 1;
    rows = (uint32_t *)malloc(cap * sizeof(uint32_t));
    cols = (uint32_t *)malloc(cap * sizeof(uint32_t));
    vals = (double *)malloc(cap * sizeof(double));
    if (!rows || !cols || !vals) {
        free(rows); free(cols); free(vals);
        return GEOMLIB_SPARSE_ERROR_MEMORY;
    }

    for (tri = 0; tri < mesh->triangle_count; ++tri) {
        uint32_t i0 = mesh->indices[tri * 3 + 0];
        uint32_t i1 = mesh->indices[tri * 3 + 1];
        uint32_t i2 = mesh->indices[tri * 3 + 2];
        double w0 = cotan_at(&mesh->vertices[i1], &mesh->vertices[i0], &mesh->vertices[i2]);
        double w1 = cotan_at(&mesh->vertices[i0], &mesh->vertices[i1], &mesh->vertices[i2]);
        double w2 = cotan_at(&mesh->vertices[i0], &mesh->vertices[i2], &mesh->vertices[i1]);
        rows[ntrip] = i0; cols[ntrip] = i1; vals[ntrip] = -0.5 * w2; ++ntrip;
        rows[ntrip] = i1; cols[ntrip] = i0; vals[ntrip] = -0.5 * w2; ++ntrip;
        rows[ntrip] = i0; cols[ntrip] = i2; vals[ntrip] = -0.5 * w1; ++ntrip;
        rows[ntrip] = i2; cols[ntrip] = i0; vals[ntrip] = -0.5 * w1; ++ntrip;
        rows[ntrip] = i1; cols[ntrip] = i2; vals[ntrip] = -0.5 * w0; ++ntrip;
        rows[ntrip] = i2; cols[ntrip] = i1; vals[ntrip] = -0.5 * w0; ++ntrip;
        rows[ntrip] = i0; cols[ntrip] = i0; vals[ntrip] = 0.5 * (w1 + w2); ++ntrip;
        rows[ntrip] = i1; cols[ntrip] = i1; vals[ntrip] = 0.5 * (w0 + w2); ++ntrip;
        rows[ntrip] = i2; cols[ntrip] = i2; vals[ntrip] = 0.5 * (w0 + w1); ++ntrip;
    }

    rc = geomlib_spmat_from_triplets(rows, cols, vals, ntrip, n, n, out);
    free(rows); free(cols); free(vals);
    return rc;
}

int geomlib_massmatrix(const geomlib_trimesh *mesh, geomlib_spmat *out)
{
    size_t tri, ntrip = 0, cap, n;
    uint32_t *rows, *cols;
    double *vals;
    int rc;

    if (!mesh || !out || !mesh->indices) return GEOMLIB_SPARSE_ERROR_INVALID;
    n = mesh->vertex_count;
    cap = mesh->triangle_count * 3 + 1;
    rows = (uint32_t *)malloc(cap * sizeof(uint32_t));
    cols = (uint32_t *)malloc(cap * sizeof(uint32_t));
    vals = (double *)malloc(cap * sizeof(double));
    if (!rows || !cols || !vals) {
        free(rows); free(cols); free(vals);
        return GEOMLIB_SPARSE_ERROR_MEMORY;
    }

    for (tri = 0; tri < mesh->triangle_count; ++tri) {
        uint32_t i0 = mesh->indices[tri * 3 + 0];
        uint32_t i1 = mesh->indices[tri * 3 + 1];
        uint32_t i2 = mesh->indices[tri * 3 + 2];
        geomlib_vec3 e1 = {
            mesh->vertices[i1].x - mesh->vertices[i0].x,
            mesh->vertices[i1].y - mesh->vertices[i0].y,
            mesh->vertices[i1].z - mesh->vertices[i0].z
        };
        geomlib_vec3 e2 = {
            mesh->vertices[i2].x - mesh->vertices[i0].x,
            mesh->vertices[i2].y - mesh->vertices[i0].y,
            mesh->vertices[i2].z - mesh->vertices[i0].z
        };
        double cx = (double)e1.y * e2.z - (double)e1.z * e2.y;
        double cy = (double)e1.z * e2.x - (double)e1.x * e2.z;
        double cz = (double)e1.x * e2.y - (double)e1.y * e2.x;
        double area = 0.5 * sqrt(cx * cx + cy * cy + cz * cz);
        rows[ntrip] = i0; cols[ntrip] = i0; vals[ntrip] = area / 3.0; ++ntrip;
        rows[ntrip] = i1; cols[ntrip] = i1; vals[ntrip] = area / 3.0; ++ntrip;
        rows[ntrip] = i2; cols[ntrip] = i2; vals[ntrip] = area / 3.0; ++ntrip;
    }

    rc = geomlib_spmat_from_triplets(rows, cols, vals, ntrip, n, n, out);
    free(rows); free(cols); free(vals);
    return rc;
}
