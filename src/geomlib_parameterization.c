/*
 * geomlib - Mesh Parameterization
 *
 * Wave 1: cotangent/Tutte harmonic Dirichlet solves, free-boundary LSCM
 * (complex least squares) and ARAP UV (local/global), all built on the
 * C11 sparse solver module.
 *
 * Papers (clean-room, no library source consulted):
 *  - Levy, Petitjean, Ray, Maillot, "Least Squares Conformal Maps for
 *    Automatic Texture Atlas Generation", SIGGRAPH 2002.
 *  - Liu, Zhang, Gotsman, Gu, "A Local/Global Approach to Mesh
 *    Parameterization", SGP 2008.
 *  - Tutte, "How to draw a graph", Proc. LMS 1963.
 */

#include "geomlib/geomlib_parameterization.h"
#include "geomlib/geomlib_healing.h"
#include "geomlib/geomlib_sparse.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------------- */
/* Boundary helpers (kept from the original API)                             */
/* ------------------------------------------------------------------------- */

typedef struct { uint32_t a, b; } boundary_edge;

static int compare_boundary_edge(const void *pa, const void *pb)
{
    const boundary_edge *x = (const boundary_edge *)pa;
    const boundary_edge *y = (const boundary_edge *)pb;
    if (x->a < y->a) return -1;
    if (x->a > y->a) return 1;
    if (x->b < y->b) return -1;
    if (x->b > y->b) return 1;
    return 0;
}

int geomlib_find_boundary_loop(
    const geomlib_trimesh *mesh,
    uint32_t **boundary_vertices,
    size_t *boundary_count
) {
    size_t nedge, i, nboundary = 0, count = 0;
    boundary_edge *edges;
    uint32_t *boundary, *degree, *adj;
    size_t *offset;
    uint32_t start, cur, prev;

    if (!mesh || !boundary_vertices || !boundary_count || !mesh->indices)
        return GEOMLIB_MESH_ERROR_INVALID;
    nedge = mesh->triangle_count * 3;
    edges = (boundary_edge *)malloc((nedge ? nedge : 1) * sizeof(boundary_edge));
    degree = (uint32_t *)calloc(mesh->vertex_count, sizeof(uint32_t));
    if (!edges || !degree) {
        free(edges); free(degree);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < nedge; ++i) {
        uint32_t a = mesh->indices[i];
        uint32_t b = mesh->indices[(i / 3) * 3 + ((i + 1) % 3)];
        if (a > b) { uint32_t t = a; a = b; b = t; }
        edges[i].a = a;
        edges[i].b = b;
    }
    qsort(edges, nedge, sizeof(boundary_edge), compare_boundary_edge);

    /* Keep edges with odd multiplicity (normally exactly one incident face). */
    for (i = 0; i < nedge;) {
        size_t j = i + 1;
        while (j < nedge && compare_boundary_edge(&edges[i], &edges[j]) == 0) ++j;
        if (((j - i) & 1u) != 0) {
            edges[nboundary++] = edges[i];
            degree[edges[i].a]++;
            degree[edges[i].b]++;
        }
        i = j;
    }
    if (nboundary == 0) {
        free(edges); free(degree);
        *boundary_vertices = NULL;
        *boundary_count = 0;
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    offset = (size_t *)calloc(mesh->vertex_count + 1, sizeof(size_t));
    adj = (uint32_t *)malloc(nboundary * 2 * sizeof(uint32_t));
    boundary = (uint32_t *)malloc(mesh->vertex_count * sizeof(uint32_t));
    if (!offset || !adj || !boundary) {
        free(edges); free(degree); free(offset); free(adj); free(boundary);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < mesh->vertex_count; ++i) offset[i+1] = offset[i] + degree[i];
    {
        size_t *cursor = (size_t *)malloc(mesh->vertex_count * sizeof(size_t));
        if (!cursor) {
            free(edges); free(degree); free(offset); free(adj); free(boundary);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        memcpy(cursor, offset, mesh->vertex_count * sizeof(size_t));
        for (i = 0; i < nboundary; ++i) {
            uint32_t a = edges[i].a, b = edges[i].b;
            adj[cursor[a]++] = b;
            adj[cursor[b]++] = a;
        }
        free(cursor);
    }

    start = edges[0].a;
    cur = start;
    prev = UINT32_MAX;
    do {
        if (count >= mesh->vertex_count) break;
        boundary[count++] = cur;
        uint32_t next = UINT32_MAX;
        for (i = offset[cur]; i < offset[cur+1]; ++i) {
            if (adj[i] != prev) { next = adj[i]; break; }
        }
        if (next == UINT32_MAX) break;
        prev = cur;
        cur = next;
    } while (cur != start);

    free(edges); free(degree); free(offset); free(adj);
    *boundary_vertices = boundary;
    *boundary_count = count;
    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_pin_boundary_circle(const uint32_t *boundary_vertices,
                                 size_t boundary_count, geomlib_vec2 *uvs)
{
    for (size_t i = 0; i < boundary_count; ++i) {
        float angle = 2.0f * (float)M_PI * (float)i / (float)boundary_count;
        uint32_t v = boundary_vertices[i];
        uvs[v].x = cosf(angle) * 0.5f + 0.5f;
        uvs[v].y = sinf(angle) * 0.5f + 0.5f;
    }
}

void geomlib_pin_boundary_square(const uint32_t *boundary_vertices,
                                 size_t boundary_count, geomlib_vec2 *uvs)
{
    for (size_t i = 0; i < boundary_count; ++i) {
        float t = (float)i / (float)boundary_count;
        uint32_t v = boundary_vertices[i];
        if (t < 0.25f) { uvs[v].x = t * 4.0f; uvs[v].y = 0.0f; }
        else if (t < 0.5f) { uvs[v].x = 1.0f; uvs[v].y = (t - 0.25f) * 4.0f; }
        else if (t < 0.75f) { uvs[v].x = 1.0f - (t - 0.5f) * 4.0f; uvs[v].y = 1.0f; }
        else { uvs[v].x = 0.0f; uvs[v].y = 1.0f - (t - 0.75f) * 4.0f; }
    }
}

/* ------------------------------------------------------------------------- */
/* Dirichlet solver on a sparse SPD matrix                                   */
/* ------------------------------------------------------------------------- */

/* Solve L_II x_I = -L_IB x_B with x_B = values on fixed vertices.
 * Fixed rows/columns are replaced by identity so the system stays SPD. */
static int solve_dirichlet(const geomlib_spmat *l, const bool *fixed,
                           const double *values, double *x)
{
    size_t n = l->nrows;
    size_t cap = l->nnz + n + 1;
    uint32_t *rows = (uint32_t *)malloc(cap * sizeof(uint32_t));
    uint32_t *cols = (uint32_t *)malloc(cap * sizeof(uint32_t));
    double *vals = (double *)malloc(cap * sizeof(double));
    double *b = (double *)calloc(n, sizeof(double));
    geomlib_spmat a;
    geomlib_sparse_solve_params params = geomlib_sparse_solve_params_default();
    size_t k, count = 0;
    int rc;

    if (!rows || !cols || !vals || !b) {
        free(rows); free(cols); free(vals); free(b);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (size_t r = 0; r < n; ++r) {
        if (fixed[r]) {
            rows[count] = (uint32_t)r;
            cols[count] = (uint32_t)r;
            vals[count++] = 1.0;
            b[r] = values[r];
            continue;
        }
        b[r] = 0.0;
        for (k = l->rowptr[r]; k < l->rowptr[r + 1]; ++k) {
            uint32_t c = l->colidx[k];
            if (fixed[c]) {
                b[r] -= l->values[k] * values[c];
            } else {
                rows[count] = (uint32_t)r;
                cols[count] = c;
                vals[count++] = l->values[k];
            }
        }
    }

    params.max_iterations = (int)(20 * n + 100);
    params.tolerance = 1e-11;
    rc = geomlib_spmat_from_triplets(rows, cols, vals, count, n, n, &a);
    if (rc == GEOMLIB_SPARSE_SUCCESS) {
        rc = geomlib_solve_cg_sparse(&a, b, x, &params);
        geomlib_spmat_free(&a);
    }
    free(rows); free(cols); free(vals); free(b);
    return rc == GEOMLIB_SPARSE_SUCCESS ? GEOMLIB_MESH_SUCCESS
                                        : GEOMLIB_MESH_ERROR_INVALID;
}

/* ------------------------------------------------------------------------- */
/* Cotangent / Tutte harmonic parameterization                               */
/* ------------------------------------------------------------------------- */

int geomlib_parameterize_harmonic(const geomlib_trimesh *mesh,
                                  geomlib_boundary_shape boundary,
                                  geomlib_uv_coords *uvs)
{
    uint32_t *boundary_verts = NULL;
    size_t boundary_count = 0;
    bool *fixed;
    geomlib_spmat l;
    double *u, *v;
    int rc;

    if (!mesh || !uvs || !uvs->uvs) return GEOMLIB_MESH_ERROR_INVALID;
    rc = geomlib_find_boundary_loop(mesh, &boundary_verts, &boundary_count);
    if (rc != GEOMLIB_MESH_SUCCESS || boundary_count == 0) return GEOMLIB_MESH_ERROR_INVALID;

    fixed = (bool *)calloc(mesh->vertex_count, sizeof(bool));
    u = (double *)calloc(mesh->vertex_count, sizeof(double));
    v = (double *)calloc(mesh->vertex_count, sizeof(double));
    if (!fixed || !u || !v) {
        free(fixed); free(u); free(v); free(boundary_verts);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t i = 0; i < boundary_count; ++i) fixed[boundary_verts[i]] = true;
    if (boundary == GEOMLIB_BOUNDARY_CIRCLE)
        geomlib_pin_boundary_circle(boundary_verts, boundary_count, uvs->uvs);
    else
        geomlib_pin_boundary_square(boundary_verts, boundary_count, uvs->uvs);
    for (size_t i = 0; i < mesh->vertex_count; ++i) {
        u[i] = uvs->uvs[i].x;
        v[i] = uvs->uvs[i].y;
    }

    rc = geomlib_cotmatrix(mesh, &l);
    if (rc == GEOMLIB_SPARSE_SUCCESS) {
        rc = solve_dirichlet(&l, fixed, u, u);
        if (rc == GEOMLIB_MESH_SUCCESS) rc = solve_dirichlet(&l, fixed, v, v);
        geomlib_spmat_free(&l);
    }
    if (rc == GEOMLIB_MESH_SUCCESS) {
        for (size_t i = 0; i < mesh->vertex_count; ++i) {
            uvs->uvs[i].x = (float)u[i];
            uvs->uvs[i].y = (float)v[i];
        }
    }

    free(fixed); free(u); free(v); free(boundary_verts);
    return rc;
}

/* ------------------------------------------------------------------------- */
/* Sparse least-squares helper (used by LSCM and ARAP)                       */
/* ------------------------------------------------------------------------- */

typedef struct {
    uint32_t *rows, *cols;
    double *vals;
    double *b;
    size_t count, cap;
    size_t ncols;
} ls_system;

static int ls_init(ls_system *s, size_t ncols, size_t hint)
{
    s->cap = hint ? hint : 256;
    s->rows = (uint32_t *)malloc(s->cap * sizeof(uint32_t));
    s->cols = (uint32_t *)malloc(s->cap * sizeof(uint32_t));
    s->vals = (double *)malloc(s->cap * sizeof(double));
    s->b = NULL;
    s->count = 0;
    s->ncols = ncols;
    if (!s->rows || !s->cols || !s->vals) {
        free(s->rows); free(s->cols); free(s->vals);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    return GEOMLIB_MESH_SUCCESS;
}

static int ls_add(ls_system *s, uint32_t row, uint32_t col, double val)
{
    if (s->count == s->cap) {
        size_t nc = s->cap * 2;
        uint32_t *r = (uint32_t *)realloc(s->rows, nc * sizeof(uint32_t));
        uint32_t *c = (uint32_t *)realloc(s->cols, nc * sizeof(uint32_t));
        double *v = (double *)realloc(s->vals, nc * sizeof(double));
        if (!r || !c || !v) { free(r); free(c); free(v); return GEOMLIB_MESH_ERROR_MEMORY; }
        s->rows = r; s->cols = c; s->vals = v; s->cap = nc;
    }
    s->rows[s->count] = row;
    s->cols[s->count] = col;
    s->vals[s->count++] = val;
    return GEOMLIB_MESH_SUCCESS;
}

static int ls_solve(ls_system *s, size_t nrows, double *x)
{
    geomlib_spmat a;
    geomlib_sparse_solve_params params = geomlib_sparse_solve_params_default();
    int rc = geomlib_spmat_from_triplets(s->rows, s->cols, s->vals, s->count,
                                         nrows, s->ncols, &a);
    if (rc != GEOMLIB_SPARSE_SUCCESS) return GEOMLIB_MESH_ERROR_INVALID;
    params.max_iterations = (int)(50 * s->ncols + 200);
    params.tolerance = 1e-10;
    rc = geomlib_solve_cgls(&a, s->b, x, &params);
    geomlib_spmat_free(&a);
    return rc == GEOMLIB_SPARSE_SUCCESS ? GEOMLIB_MESH_SUCCESS
                                        : GEOMLIB_MESH_ERROR_INVALID;
}

static void ls_free(ls_system *s)
{
    free(s->rows); free(s->cols); free(s->vals); free(s->b);
    memset(s, 0, sizeof(*s));
}

/* ------------------------------------------------------------------------- */
/* Free-boundary LSCM                                                        */
/* ------------------------------------------------------------------------- */

int geomlib_parameterize_lscm_free(const geomlib_trimesh *mesh,
                                   uint32_t pin_a, uint32_t pin_b,
                                   geomlib_uv_coords *uvs)
{
    uint32_t *bverts = NULL;
    size_t bcount = 0;
    ls_system s;
    size_t nrows, t, i;
    double *x;
    int rc;

    if (!mesh || !uvs || !uvs->uvs || mesh->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    if (pin_a == UINT32_MAX || pin_b == UINT32_MAX ||
        pin_a >= mesh->vertex_count || pin_b >= mesh->vertex_count ||
        pin_a == pin_b) {
        rc = geomlib_find_boundary_loop(mesh, &bverts, &bcount);
        if (rc != GEOMLIB_MESH_SUCCESS || bcount < 2) {
            free(bverts);
            return GEOMLIB_MESH_ERROR_INVALID;
        }
        /* Pick the farthest-apart boundary pair for stability. */
        {
            float best = -1.0f;
            size_t p = 0, q = 0;
            for (size_t a = 0; a < bcount; ++a) {
                for (size_t b = a + 1; b < bcount; ++b) {
                    geomlib_vec3 d = geomlib_vec3_sub(mesh->vertices[bverts[a]],
                                                      mesh->vertices[bverts[b]]);
                    float dist = geomlib_vec3_length(d);
                    if (dist > best) { best = dist; p = a; q = b; }
                }
            }
            pin_a = bverts[p];
            pin_b = bverts[q];
        }
    }

    nrows = mesh->triangle_count * 2 + 4;
    if (ls_init(&s, (uint32_t)(mesh->vertex_count * 2),
                mesh->triangle_count * 12 + 8) != GEOMLIB_MESH_SUCCESS) {
        free(bverts);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    s.b = (double *)calloc(nrows, sizeof(double));
    x = (double *)calloc(mesh->vertex_count * 2, sizeof(double));
    if (!s.b || !x) {
        ls_free(&s); free(x); free(bverts);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (t = 0; t < mesh->triangle_count; ++t) {
        uint32_t v[3];
        geomlib_vec3 p0, e1, e2, n, t1, t2;
        double z1r, z1i, z2r, z2i, denom;
        double c0r, c0i, c1r, c1i, c2r, c2i;
        uint32_t re = (uint32_t)(2 * t), im = re + 1;

        v[0] = mesh->indices[t * 3 + 0];
        v[1] = mesh->indices[t * 3 + 1];
        v[2] = mesh->indices[t * 3 + 2];
        p0 = mesh->vertices[v[0]];
        e1 = geomlib_vec3_sub(mesh->vertices[v[1]], p0);
        e2 = geomlib_vec3_sub(mesh->vertices[v[2]], p0);
        n = geomlib_vec3_cross(e1, e2);
        if (geomlib_vec3_length(n) < 1e-12f) continue;
        n = geomlib_vec3_normalize(n);
        t1 = geomlib_vec3_normalize(e1);
        t2 = geomlib_vec3_cross(n, t1);
        z1r = geomlib_vec3_dot(e1, t1);
        z1i = geomlib_vec3_dot(e1, t2);
        z2r = geomlib_vec3_dot(e2, t1);
        z2i = geomlib_vec3_dot(e2, t2);
        denom = sqrt((z1r*z1r + z1i*z1i + z2r*z2r + z2i*z2i) * 0.5);
        if (denom < 1e-12) continue;
        z1r /= denom; z1i /= denom; z2r /= denom; z2i /= denom;

        /* (z1 - z2) w0 + z2 w1 - z1 w2 = 0 */
        c0r = z1r - z2r; c0i = z1i - z2i;
        c1r = z2r;       c1i = z2i;
        c2r = -z1r;      c2i = -z1i;
        {
            double cr[3] = {c0r, c1r, c2r};
            double ci[3] = {c0i, c1i, c2i};
            for (i = 0; i < 3; ++i) {
                if (ls_add(&s, re, 2*v[i], cr[i]) != GEOMLIB_MESH_SUCCESS ||
                    ls_add(&s, re, 2*v[i]+1, -ci[i]) != GEOMLIB_MESH_SUCCESS ||
                    ls_add(&s, im, 2*v[i], ci[i]) != GEOMLIB_MESH_SUCCESS ||
                    ls_add(&s, im, 2*v[i]+1, cr[i]) != GEOMLIB_MESH_SUCCESS) {
                    ls_free(&s); free(x); free(bverts);
                    return GEOMLIB_MESH_ERROR_MEMORY;
                }
            }
        }
    }

    /* Pin two vertices: w_a = 0, w_b = 1. Heavy rows make them exact. */
    {
        uint32_t row = (uint32_t)(mesh->triangle_count * 2);
        const double weight = 1000.0;
        ls_add(&s, row, 2*pin_a, weight); s.b[row] = 0.0;
        ls_add(&s, row+1, 2*pin_a+1, weight); s.b[row+1] = 0.0;
        ls_add(&s, row+2, 2*pin_b, weight); s.b[row+2] = weight;
        ls_add(&s, row+3, 2*pin_b+1, weight); s.b[row+3] = 0.0;
    }

    rc = ls_solve(&s, nrows, x);
    if (rc == GEOMLIB_MESH_SUCCESS) {
        for (i = 0; i < mesh->vertex_count; ++i) {
            uvs->uvs[i].x = (float)x[2*i];
            uvs->uvs[i].y = (float)x[2*i+1];
        }
    }
    ls_free(&s); free(x); free(bverts);
    return rc;
}

/* ------------------------------------------------------------------------- */
/* ARAP parameterization (local/global)                                      */
/* ------------------------------------------------------------------------- */

typedef struct {
    size_t count, cap;
    uint32_t *neighbor;
    double *ax, *ay;     /* local 3D coordinates relative to the center */
} one_ring;

static void ring_free(one_ring *r)
{
    free(r->neighbor); free(r->ax); free(r->ay);
    memset(r, 0, sizeof(*r));
}

static int build_rings(const geomlib_trimesh *mesh, one_ring *rings)
{
    size_t v, t, i;
    memset(rings, 0, mesh->vertex_count * sizeof(one_ring));
    /* Count incidences first. */
    size_t *counts = (size_t *)calloc(mesh->vertex_count, sizeof(size_t));
    if (!counts) return GEOMLIB_MESH_ERROR_MEMORY;
    for (t = 0; t < mesh->triangle_count; ++t) {
        for (i = 0; i < 3; ++i) counts[mesh->indices[t*3+i]] += 2;
    }
    for (v = 0; v < mesh->vertex_count; ++v) {
        rings[v].cap = counts[v] + 1;
        rings[v].neighbor = (uint32_t *)malloc(rings[v].cap * sizeof(uint32_t));
        rings[v].ax = (double *)malloc(rings[v].cap * sizeof(double));
        rings[v].ay = (double *)malloc(rings[v].cap * sizeof(double));
        if (!rings[v].neighbor || !rings[v].ax || !rings[v].ay) {
            for (i = 0; i <= v; ++i) ring_free(&rings[i]);
            free(counts);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        rings[v].count = 0;
    }
    free(counts);

    /* Compute ONE stable tangent frame per vertex. Recomputing it per
     * incident triangle would make the local rest coordinates inconsistent
     * and systematically drift the ARAP solution. */
    {
        geomlib_vec3 *normals = (geomlib_vec3 *)calloc(mesh->vertex_count, sizeof(geomlib_vec3));
        geomlib_vec3 *t1s = (geomlib_vec3 *)calloc(mesh->vertex_count, sizeof(geomlib_vec3));
        geomlib_vec3 *t2s = (geomlib_vec3 *)calloc(mesh->vertex_count, sizeof(geomlib_vec3));
        if (!normals || !t1s || !t2s) {
            free(normals); free(t1s); free(t2s);
            for (v = 0; v < mesh->vertex_count; ++v) ring_free(&rings[v]);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        for (t = 0; t < mesh->triangle_count; ++t) {
            uint32_t a = mesh->indices[t*3], b = mesh->indices[t*3+1], c = mesh->indices[t*3+2];
            geomlib_vec3 n = geomlib_vec3_cross(
                geomlib_vec3_sub(mesh->vertices[b], mesh->vertices[a]),
                geomlib_vec3_sub(mesh->vertices[c], mesh->vertices[a]));
            normals[a] = geomlib_vec3_add(normals[a], n);
            normals[b] = geomlib_vec3_add(normals[b], n);
            normals[c] = geomlib_vec3_add(normals[c], n);
        }
        for (v = 0; v < mesh->vertex_count; ++v) {
            geomlib_vec3 nrm = geomlib_vec3_normalize(normals[v]);
            geomlib_vec3 ref = geomlib_vec3_create(1, 0, 0);
            if (geomlib_vec3_length(nrm) < 0.5f) nrm = geomlib_vec3_create(0, 0, 1);
            t1s[v] = geomlib_vec3_sub(ref, geomlib_vec3_mul(nrm, geomlib_vec3_dot(ref, nrm)));
            if (geomlib_vec3_length(t1s[v]) < 1e-6f) {
                ref = geomlib_vec3_create(0, 1, 0);
                t1s[v] = geomlib_vec3_sub(ref, geomlib_vec3_mul(nrm, geomlib_vec3_dot(ref, nrm)));
            }
            t1s[v] = geomlib_vec3_normalize(t1s[v]);
            t2s[v] = geomlib_vec3_cross(nrm, t1s[v]);
        }
        for (t = 0; t < mesh->triangle_count; ++t) {
            for (i = 0; i < 3; ++i) {
                uint32_t center = mesh->indices[t*3+i];
                uint32_t next = mesh->indices[t*3 + (i+1)%3];
                uint32_t prev = mesh->indices[t*3 + (i+2)%3];
                int k;
                for (k = 0; k < 2; ++k) {
                    uint32_t nb = (k == 0) ? next : prev;
                    geomlib_vec3 dv = geomlib_vec3_sub(mesh->vertices[nb], mesh->vertices[center]);
                    one_ring *r = &rings[center];
                    if (r->count >= r->cap) break;
                    r->neighbor[r->count] = nb;
                    /* The boundary is pinned into [0,1]^2, so rest edges use
                     * half the 3D length to match that scale. */
                    r->ax[r->count] = 0.5 * geomlib_vec3_dot(dv, t1s[center]);
                    r->ay[r->count] = 0.5 * geomlib_vec3_dot(dv, t2s[center]);
                    r->count++;
                }
            }
        }
        free(normals); free(t1s); free(t2s);
    }
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_parameterize_arap(const geomlib_trimesh *mesh,
                              geomlib_boundary_shape boundary,
                              int max_iterations, geomlib_uv_coords *uvs)
{
    uint32_t *bverts = NULL;
    size_t bcount = 0;
    bool *fixed;
    one_ring *rings;
    double *cos_r, *sin_r;
    double *target_u, *target_v;
    int rc, iter;
    size_t v, i;

    if (!mesh || !uvs || !uvs->uvs || mesh->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    rc = geomlib_find_boundary_loop(mesh, &bverts, &bcount);
    if (rc != GEOMLIB_MESH_SUCCESS || bcount < 3) {
        free(bverts);
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    if (max_iterations <= 0) max_iterations = 10;

    fixed = (bool *)calloc(mesh->vertex_count, sizeof(bool));
    rings = (one_ring *)calloc(mesh->vertex_count, sizeof(one_ring));
    cos_r = (double *)calloc(mesh->vertex_count, sizeof(double));
    sin_r = (double *)calloc(mesh->vertex_count, sizeof(double));
    target_u = (double *)calloc(mesh->vertex_count, sizeof(double));
    target_v = (double *)calloc(mesh->vertex_count, sizeof(double));
    if (!fixed || !rings || !cos_r || !sin_r || !target_u || !target_v) {
        free(fixed); free(rings); free(cos_r); free(sin_r);
        free(target_u); free(target_v); free(bverts);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < bcount; ++i) fixed[bverts[i]] = true;
    if (boundary == GEOMLIB_BOUNDARY_CIRCLE)
        geomlib_pin_boundary_circle(bverts, bcount, uvs->uvs);
    else
        geomlib_pin_boundary_square(bverts, bcount, uvs->uvs);
    for (v = 0; v < mesh->vertex_count; ++v) {
        target_u[v] = uvs->uvs[v].x;
        target_v[v] = uvs->uvs[v].y;
    }

    /* Initialize with the cotangent harmonic solution. */
    rc = geomlib_parameterize_harmonic(mesh, boundary, uvs);
    if (rc != GEOMLIB_MESH_SUCCESS) goto arap_done;
    rc = build_rings(mesh, rings);
    if (rc != GEOMLIB_MESH_SUCCESS) goto arap_done;

    for (iter = 0; iter < max_iterations; ++iter) {
        ls_system s;
        double *x = (double *)calloc(mesh->vertex_count * 2, sizeof(double));
        size_t nrows;
        if (!x) { rc = GEOMLIB_MESH_ERROR_MEMORY; break; }

        /* Local step: best-fit 2D rotation per vertex (complex unit phase). */
        for (v = 0; v < mesh->vertex_count; ++v) {
            one_ring *r = &rings[v];
            double re = 0.0, im = 0.0, len;
            for (i = 0; i < r->count; ++i) {
                uint32_t nb = r->neighbor[i];
                double bu = uvs->uvs[nb].x - uvs->uvs[v].x;
                double bv = uvs->uvs[nb].y - uvs->uvs[v].y;
                double au = r->ax[i], av = r->ay[i];
                re += bu * au + bv * av;
                im += bv * au - bu * av;
            }
            len = sqrt(re*re + im*im);
            if (len > 1e-12) { cos_r[v] = re / len; sin_r[v] = im / len; }
            else { cos_r[v] = 1.0; sin_r[v] = 0.0; }
        }

        /* Global step: solve for w minimizing the ARAP residual. */
        nrows = 0;
        for (v = 0; v < mesh->vertex_count; ++v) nrows += 2 * rings[v].count;
        nrows += 2 * bcount;
        rc = ls_init(&s, (uint32_t)(mesh->vertex_count * 2), nrows * 2 + 8);
        if (rc != GEOMLIB_MESH_SUCCESS) { free(x); break; }
        s.b = (double *)calloc(nrows, sizeof(double));
        if (!s.b) { ls_free(&s); free(x); rc = GEOMLIB_MESH_ERROR_MEMORY; break; }
        {
            size_t row = 0;
            for (v = 0; v < mesh->vertex_count; ++v) {
                one_ring *r = &rings[v];
                for (i = 0; i < r->count; ++i) {
                    uint32_t nb = r->neighbor[i];
                    double au = r->ax[i], av = r->ay[i];
                    double ru = cos_r[v]*au - sin_r[v]*av;
                    double rv = sin_r[v]*au + cos_r[v]*av;
                    ls_add(&s, (uint32_t)row, 2*nb, 1.0);
                    ls_add(&s, (uint32_t)row, 2*v, -1.0);
                    s.b[row] = ru; row++;
                    ls_add(&s, (uint32_t)row, 2*nb+1, 1.0);
                    ls_add(&s, (uint32_t)row, 2*v+1, -1.0);
                    s.b[row] = rv; row++;
                }
            }
            for (i = 0; i < bcount; ++i) {
                uint32_t bv = bverts[i];
                /* Boundary constraints must dominate the many edge residuals;
                 * a weak pin lets the whole map drift outside the disk. */
                const double weight = 100.0;
                ls_add(&s, (uint32_t)row, 2*bv, weight);
                s.b[row] = weight * target_u[bv]; row++;
                ls_add(&s, (uint32_t)row, 2*bv+1, weight);
                s.b[row] = weight * target_v[bv]; row++;
            }
        }
        rc = ls_solve(&s, nrows, x);
        if (rc == GEOMLIB_MESH_SUCCESS) {
            for (v = 0; v < mesh->vertex_count; ++v) {
                uvs->uvs[v].x = (float)x[2*v];
                uvs->uvs[v].y = (float)x[2*v+1];
            }
        }
        ls_free(&s);
        free(x);
        if (rc != GEOMLIB_MESH_SUCCESS) break;
    }

arap_done:
    for (v = 0; v < mesh->vertex_count; ++v) ring_free(&rings[v]);
    free(rings); free(fixed); free(cos_r); free(sin_r);
    free(target_u); free(target_v); free(bverts);
    return rc;
}

/* ------------------------------------------------------------------------- */
/* Dispatcher and UV management                                              */
/* ------------------------------------------------------------------------- */

/* Public LSCM entry point. A pinned boundary uses the cotangent/Tutte solve;
 * a natural (free) boundary uses the two-pin complex least-squares solver. */
int geomlib_parameterize_lscm(const geomlib_trimesh *mesh,
                              geomlib_boundary_shape boundary,
                              geomlib_uv_coords *uvs)
{
    if (!mesh || !uvs) return GEOMLIB_MESH_ERROR_INVALID;
    if (boundary == GEOMLIB_BOUNDARY_NATURAL)
        return geomlib_parameterize_lscm_free(mesh, UINT32_MAX, UINT32_MAX, uvs);
    return geomlib_parameterize_harmonic(mesh, boundary, uvs);
}

geomlib_param_params geomlib_param_params_default(geomlib_param_method method)
{
    geomlib_param_params p;
    p.method = method;
    p.boundary = GEOMLIB_BOUNDARY_CIRCLE;
    p.pin_vertex_0 = UINT32_MAX;
    p.pin_vertex_1 = UINT32_MAX;
    p.max_iterations = 10;
    p.tolerance = 1e-6f;
    p.normalize_uvs = true;
    return p;
}

int geomlib_parameterize(const geomlib_trimesh *mesh,
                         const geomlib_param_params *params,
                         geomlib_uv_coords *uvs)
{
    int rc;
    if (!mesh || !params || !uvs) return GEOMLIB_MESH_ERROR_INVALID;
    switch (params->method) {
        case GEOMLIB_PARAM_LSCM:
            if (params->boundary == GEOMLIB_BOUNDARY_NATURAL)
                rc = geomlib_parameterize_lscm_free(mesh, params->pin_vertex_0,
                                                    params->pin_vertex_1, uvs);
            else
                rc = geomlib_parameterize_harmonic(mesh, params->boundary, uvs);
            break;
        case GEOMLIB_PARAM_HARMONIC:
            rc = geomlib_parameterize_harmonic(mesh, params->boundary, uvs);
            break;
        case GEOMLIB_PARAM_ARAP:
            rc = geomlib_parameterize_arap(mesh, params->boundary,
                                           params->max_iterations, uvs);
            break;
        default:
            return GEOMLIB_MESH_ERROR_INVALID;
    }
    if (rc == GEOMLIB_MESH_SUCCESS && params->normalize_uvs)
        geomlib_normalize_uvs(uvs);
    return rc;
}

void geomlib_normalize_uvs(geomlib_uv_coords *uvs)
{
    if (!uvs || !uvs->uvs || uvs->vertex_count == 0) return;
    float min_u = FLT_MAX, min_v = FLT_MAX, max_u = -FLT_MAX, max_v = -FLT_MAX;
    for (size_t i = 0; i < uvs->vertex_count; ++i) {
        if (uvs->uvs[i].x < min_u) min_u = uvs->uvs[i].x;
        if (uvs->uvs[i].y < min_v) min_v = uvs->uvs[i].y;
        if (uvs->uvs[i].x > max_u) max_u = uvs->uvs[i].x;
        if (uvs->uvs[i].y > max_v) max_v = uvs->uvs[i].y;
    }
    float su = (max_u > min_u) ? 1.0f / (max_u - min_u) : 1.0f;
    float sv = (max_v > min_v) ? 1.0f / (max_v - min_v) : 1.0f;
    for (size_t i = 0; i < uvs->vertex_count; ++i) {
        uvs->uvs[i].x = (uvs->uvs[i].x - min_u) * su;
        uvs->uvs[i].y = (uvs->uvs[i].y - min_v) * sv;
    }
}

geomlib_uv_coords geomlib_uv_coords_create(size_t vertex_count, size_t triangle_count)
{
    geomlib_uv_coords u;
    u.uvs = (geomlib_vec2 *)calloc(vertex_count, sizeof(geomlib_vec2));
    u.vertex_count = vertex_count;
    u.distortion = triangle_count ? (float *)calloc(triangle_count, sizeof(float)) : NULL;
    u.triangle_count = triangle_count;
    return u;
}

void geomlib_uv_coords_free(geomlib_uv_coords *u)
{
    if (!u) return;
    free(u->uvs); free(u->distortion);
    u->uvs = NULL; u->distortion = NULL;
    u->vertex_count = u->triangle_count = 0;
}

/* ------------------------------------------------------------------------- */
/* Distortion statistics                                                     */
/* ------------------------------------------------------------------------- */

void geomlib_compute_param_distortion(const geomlib_trimesh *mesh,
                                      const geomlib_uv_coords *uvs)
{
    if (!mesh || !uvs || !uvs->distortion) return;
    for (size_t t = 0; t < mesh->triangle_count; ++t) {
        uint32_t i0 = mesh->indices[t*3], i1 = mesh->indices[t*3+1], i2 = mesh->indices[t*3+2];
        geomlib_vec3 p0 = mesh->vertices[i0], p1 = mesh->vertices[i1], p2 = mesh->vertices[i2];
        geomlib_vec2 w0 = uvs->uvs[i0], w1 = uvs->uvs[i1], w2 = uvs->uvs[i2];
        geomlib_vec3 e1 = geomlib_vec3_sub(p1, p0), e2 = geomlib_vec3_sub(p2, p0);
        geomlib_vec2 f1 = geomlib_vec2_sub(w1, w0), f2 = geomlib_vec2_sub(w2, w0);
        float a3 = geomlib_vec3_length(geomlib_vec3_cross(e1, e2)) * 0.5f;
        float a2 = fabsf(f1.x*f2.y - f1.y*f2.x) * 0.5f;
        uvs->distortion[t] = (a3 > 1e-10f) ? (a2 / a3) : 0.0f;
    }
}

void geomlib_compute_param_stats(const geomlib_trimesh *mesh,
                                 const geomlib_uv_coords *uvs,
                                 geomlib_param_stats *stats)
{
    if (!mesh || !uvs || !stats) return;
    memset(stats, 0, sizeof(*stats));
    if (mesh->triangle_count == 0) return;

    double angle_err = 0.0;
    double stretch2 = 0.0, stretch_inf = 0.0;
    double area_err = 0.0;
    size_t angles = 0;

    for (size_t t = 0; t < mesh->triangle_count; ++t) {
        uint32_t idx[3] = { mesh->indices[t*3], mesh->indices[t*3+1], mesh->indices[t*3+2] };
        geomlib_vec3 p[3];
        geomlib_vec2 w[3];
        double s3[3], s2[3];
        for (int k = 0; k < 3; ++k) {
            p[k] = mesh->vertices[idx[k]];
            w[k] = uvs->uvs[idx[k]];
        }
        for (int k = 0; k < 3; ++k) {
            geomlib_vec3 d3 = geomlib_vec3_sub(p[(k+1)%3], p[k]);
            geomlib_vec2 d2 = geomlib_vec2_sub(w[(k+1)%3], w[k]);
            s3[k] = geomlib_vec3_length(d3);
            s2[k] = sqrt(d2.x*d2.x + d2.y*d2.y);
        }
        double signed_area = (w[1].x-w[0].x)*(w[2].y-w[0].y) - (w[1].y-w[0].y)*(w[2].x-w[0].x);
        if (signed_area < 0.0) { stats->has_flips = true; stats->flipped_count++; }

        /* Per-corner angle distortion (radians), matched by sorted edge order. */
        for (int k = 0; k < 3; ++k) {
            double a3, a2;
            int j = (k+1)%3, l = (k+2)%3;
            if (s3[j] < 1e-12 || s3[l] < 1e-12 || s2[j] < 1e-12 || s2[l] < 1e-12) continue;
            /* Use the vector dot product directly for robustness. */
            geomlib_vec3 vj = geomlib_vec3_sub(p[j], p[k]);
            geomlib_vec3 vl = geomlib_vec3_sub(p[l], p[k]);
            geomlib_vec2 uj = geomlib_vec2_sub(w[j], w[k]);
            geomlib_vec2 ul = geomlib_vec2_sub(w[l], w[k]);
            double d3d = geomlib_vec3_dot(vj, vl) / (s3[j]*s3[l]);
            double d2d = (uj.x*ul.x + uj.y*ul.y) / (s2[j]*s2[l]);
            if (d3d < -1.0) d3d = -1.0;
            if (d3d > 1.0) d3d = 1.0;
            if (d2d < -1.0) d2d = -1.0;
            if (d2d > 1.0) d2d = 1.0;
            a3 = acos(d3d); a2 = acos(d2d);
            angle_err += (a2 - a3) * (a2 - a3);
            angles++;
        }

        /* Stretch: singular-value-like ratios from the two edge directions.
         * Using edge-length ratios gives a conservative, well-defined metric. */
        for (int k = 0; k < 3; ++k) {
            if (s3[k] < 1e-12) continue;
            double r = s2[k] / s3[k];
            stretch2 += r*r;
            if (r > stretch_inf) stretch_inf = r;
        }
        if (uvs->distortion) {
            double r = uvs->distortion[t];
            area_err += (r - 1.0)*(r - 1.0);
        }
    }

    if (angles) stats->angle_distortion = (float)sqrt(angle_err / angles);
    stats->stretch_L2 = (float)sqrt(stretch2 / (3.0 * mesh->triangle_count));
    stats->stretch_Linf = (float)stretch_inf;
    stats->area_distortion = (float)sqrt(area_err / mesh->triangle_count);
}
