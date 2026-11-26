/*
 * geomlib - Variational Shape Approximation (see geomlib_vsa.h)
 *
 * Clean-room implementation of the VSA idea from Cohen-Sfady et al.,
 * "Variational Shape Approximation", SIGGRAPH 2004: cluster triangles around
 * planar proxies minimizing squared distance, then rebuild a coarse mesh from
 * the region boundaries.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_vsa.h"
#include "geomlib/geomlib_sparse.h"
#include "geomlib/geomlib_repair.h"

#include <float.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_vsa_result geomlib_vsa_result_create(void)
{
    geomlib_vsa_result r;
    memset(&r, 0, sizeof(r));
    return r;
}

void geomlib_vsa_result_free(geomlib_vsa_result *r)
{
    if (!r) return;
    free(r->region_of_triangle);
    free(r->proxy_points);
    free(r->proxy_normals);
    memset(r, 0, sizeof(*r));
}

typedef struct { uint32_t a, b; } ekey;

static int cmp_ekey(const void *pa, const void *pb)
{
    const ekey *x = (const ekey *)pa, *y = (const ekey *)pb;
    if (x->a < y->a) return -1;
    if (x->a > y->a) return 1;
    if (x->b < y->b) return -1;
    if (x->b > y->b) return 1;
    return 0;
}

static geomlib_vec3 triangle_center(const geomlib_trimesh *m, size_t t)
{
    geomlib_vec3 a = m->vertices[m->indices[t*3]];
    geomlib_vec3 b = m->vertices[m->indices[t*3+1]];
    geomlib_vec3 c = m->vertices[m->indices[t*3+2]];
    return geomlib_vec3_mul(geomlib_vec3_add(geomlib_vec3_add(a,b),c), 1.0f/3.0f);
}

/* Fit a plane to all vertices of the assigned triangles (PCA). */
static void fit_plane(const geomlib_trimesh *m, const uint32_t *assign,
                      size_t T, uint32_t region, geomlib_vec3 *point,
                      geomlib_vec3 *normal)
{
    double c[3] = {0,0,0};
    double cov[9] = {0};
    size_t count = 0, t, k;
    for (t = 0; t < T; ++t) {
        if (assign[t] != region) continue;
        for (k = 0; k < 3; ++k) {
            geomlib_vec3 p = m->vertices[m->indices[t*3+k]];
            c[0]+=p.x; c[1]+=p.y; c[2]+=p.z;
            count++;
        }
    }
    if (!count) { *point = geomlib_vec3_create(0,0,0); *normal = geomlib_vec3_create(0,0,1); return; }
    c[0]/=count; c[1]/=count; c[2]/=count;
    for (t = 0; t < T; ++t) {
        if (assign[t] != region) continue;
        for (k = 0; k < 3; ++k) {
            geomlib_vec3 p = m->vertices[m->indices[t*3+k]];
            double v[3] = {p.x-c[0], p.y-c[1], p.z-c[2]};
            int i, j;
            for (i = 0; i < 3; ++i)
                for (j = 0; j < 3; ++j) cov[i*3+j] += v[i]*v[j];
        }
    }
    {
        double eig[3], vec[9];
        if (geomlib_eigen_sym_dense(cov, 3, eig, vec) == GEOMLIB_SPARSE_SUCCESS)
            *normal = geomlib_vec3_create((float)vec[0], (float)vec[1], (float)vec[2]);
        else
            *normal = geomlib_vec3_create(0,0,1);
    }
    *point = geomlib_vec3_create((float)c[0], (float)c[1], (float)c[2]);
}

static double plane_cost(const geomlib_trimesh *m, size_t t,
                         geomlib_vec3 p, geomlib_vec3 n)
{
    double cost = 0.0;
    int k;
    for (k = 0; k < 3; ++k) {
        geomlib_vec3 v = geomlib_vec3_sub(m->vertices[m->indices[t*3+k]], p);
        float d = geomlib_vec3_dot(v, n);
        cost += (double)d*d;
    }
    return cost;
}

typedef struct {
    geomlib_vec3 *v;
    size_t count, cap;
    uint32_t *idx;
    size_t icount, icap;
} soup;

static void add_vertex(soup *s, geomlib_vec3 p, uint32_t *out)
{
    if (s->count == s->cap) {
        size_t nc = s->cap ? s->cap*2 : 4096;
        geomlib_vec3 *nv = (geomlib_vec3 *)realloc(s->v, nc*sizeof(geomlib_vec3));
        if (!nv) return;
        s->v = nv; s->cap = nc;
    }
    s->v[s->count] = p;
    *out = (uint32_t)s->count++;
}

static void add_tri(soup *s, uint32_t a, uint32_t b, uint32_t c)
{
    if (s->icount + 3 > s->icap) {
        size_t nc = s->icap ? s->icap*2 : 4096;
        uint32_t *ni = (uint32_t *)realloc(s->idx, nc*sizeof(uint32_t));
        if (!ni) return;
        s->idx = ni; s->icap = nc;
    }
    s->idx[s->icount++] = a;
    s->idx[s->icount++] = b;
    s->idx[s->icount++] = c;
}

int geomlib_vsa_simplify(const geomlib_trimesh *input,
                         size_t region_count,
                         int iterations,
                         geomlib_vsa_result *result,
                         geomlib_trimesh *output)
{
    size_t T, V, r, t, k, i;
    uint32_t *assign, *seeds;
    geomlib_vec3 *centers, *points, *normals;
    soup s;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!input || !output || !input->vertices || !input->indices ||
        input->triangle_count == 0 || !result)
        return GEOMLIB_MESH_ERROR_INVALID;
    T = input->triangle_count;
    V = input->vertex_count;
    if (region_count < 1) region_count = 1;
    if (region_count > T) region_count = T;
    if (iterations <= 0) iterations = 10;
    if (iterations > 50) iterations = 50;

    centers = (geomlib_vec3 *)malloc(T*sizeof(geomlib_vec3));
    assign = (uint32_t *)malloc(T*sizeof(uint32_t));
    seeds = (uint32_t *)malloc(region_count*sizeof(uint32_t));
    points = (geomlib_vec3 *)malloc(region_count*sizeof(geomlib_vec3));
    normals = (geomlib_vec3 *)malloc(region_count*sizeof(geomlib_vec3));
    if (!centers || !assign || !seeds || !points || !normals) {
        free(centers); free(assign); free(seeds); free(points); free(normals);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (t = 0; t < T; ++t) centers[t] = triangle_center(input, t);

    /* Farthest-point seeding. */
    seeds[0] = 0;
    for (r = 1; r < region_count; ++r) {
        double best = -1.0;
        for (t = 0; t < T; ++t) {
            double mind = DBL_MAX;
            for (i = 0; i < r; ++i) {
                geomlib_vec3 d = geomlib_vec3_sub(centers[t], centers[seeds[i]]);
                double dd = geomlib_vec3_dot(d, d);
                if (dd < mind) mind = dd;
            }
            if (mind > best) { best = mind; seeds[r] = (uint32_t)t; }
        }
    }
    for (r = 0; r < region_count; ++r) {
        points[r] = centers[seeds[r]];
        normals[r] = geomlib_trimesh_triangle_normal(input, seeds[r]);
    }
    for (t = 0; t < T; ++t) assign[t] = 0;

    for (i = 0; i < (size_t)iterations; ++i) {
        int changed = 0;
        for (t = 0; t < T; ++t) {
            double best = DBL_MAX;
            uint32_t best_r = 0;
            for (r = 0; r < region_count; ++r) {
                double cost = plane_cost(input, t, points[r], normals[r]);
                if (cost < best) { best = cost; best_r = (uint32_t)r; }
            }
            if (assign[t] != best_r) { assign[t] = best_r; changed = 1; }
        }
        for (r = 0; r < region_count; ++r)
            fit_plane(input, assign, T, (uint32_t)r, &points[r], &normals[r]);
        if (!changed) break;
    }

    /* Fill the result clustering. */
    {
        geomlib_vsa_result res = geomlib_vsa_result_create();
        res.region_of_triangle = (uint32_t *)malloc(T*sizeof(uint32_t));
        res.proxy_points = (geomlib_vec3 *)malloc(region_count*sizeof(geomlib_vec3));
        res.proxy_normals = (geomlib_vec3 *)malloc(region_count*sizeof(geomlib_vec3));
        if (!res.region_of_triangle || !res.proxy_points || !res.proxy_normals) {
            geomlib_vsa_result_free(&res);
            free(centers); free(assign); free(seeds); free(points); free(normals);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        memcpy(res.region_of_triangle, assign, T*sizeof(uint32_t));
        memcpy(res.proxy_points, points, region_count*sizeof(geomlib_vec3));
        memcpy(res.proxy_normals, normals, region_count*sizeof(geomlib_vec3));
        res.triangle_count = T;
        res.region_count = region_count;
        *result = res;
    }

    /* Build the coarse mesh from each region's boundary edges. */
    memset(&s, 0, sizeof(s));
    for (r = 0; r < region_count; ++r) {
        uint32_t centroid_idx = 0;
        ekey *edges = (ekey *)malloc(T*3*sizeof(ekey));
        size_t ne = 0;
        if (!edges) continue;
        for (t = 0; t < T; ++t) {
            if (assign[t] != r) continue;
            for (k = 0; k < 3; ++k) {
                uint32_t a = input->indices[t*3+k], b = input->indices[t*3+(k+1)%3];
                edges[ne].a = a < b ? a : b;
                edges[ne].b = a < b ? b : a;
                ne++;
            }
        }
        qsort(edges, ne, sizeof(ekey), cmp_ekey);
        {
            size_t n = 0, j;
            for (i = 0; i < ne;) {
                j = i+1;
                while (j < ne && cmp_ekey(&edges[i], &edges[j]) == 0) ++j;
                if ((j-i) & 1u) edges[n++] = edges[i];
                i = j;
            }
            ne = n;
        }
        if (!ne) { free(edges); continue; }

        add_vertex(&s, points[r], &centroid_idx);
        /* Adjacency of boundary vertices. */
        {
            size_t *count = (size_t *)calloc(V, sizeof(size_t));
            size_t *off = (size_t *)malloc((V+1)*sizeof(size_t));
            uint32_t *adj = (uint32_t *)malloc(ne*2*sizeof(uint32_t));
            size_t *cursor;
            if (!count || !off || !adj) {
                free(count); free(off); free(adj); free(edges); continue;
            }
            for (i = 0; i < ne; ++i) { count[edges[i].a]++; count[edges[i].b]++; }
            off[0] = 0;
            for (i = 0; i < V; ++i) off[i+1] = off[i] + count[i];
            cursor = (size_t *)malloc(V*sizeof(size_t));
            if (!cursor) { free(count); free(off); free(adj); free(edges); continue; }
            memcpy(cursor, off, V*sizeof(size_t));
            for (i = 0; i < ne; ++i) {
                adj[cursor[edges[i].a]++] = edges[i].b;
                adj[cursor[edges[i].b]++] = edges[i].a;
            }
            free(cursor);
            /* Walk cycles; emit a fan from the proxy centroid. */
            {
                bool *used = (bool *)calloc(ne, sizeof(bool));
                if (used) {
                    for (i = 0; i < ne; ++i) {
                        uint32_t start, cur, prev;
                        size_t guard = 0;
                        if (used[i]) continue;
                        used[i] = true;
                        start = edges[i].a;
                        cur = edges[i].b;
                        prev = edges[i].a;
                        {
                            uint32_t va = 0, vb = 0;
                            add_vertex(&s, input->vertices[start], &va);
                            add_vertex(&s, input->vertices[cur], &vb);
                            add_tri(&s, centroid_idx, va, vb);
                        }
                        while (cur != start && guard++ < ne+2) {
                            uint32_t next = UINT32_MAX;
                            size_t j;
                            for (j = off[cur]; j < off[cur+1]; ++j) {
                                if (adj[j] != prev) { next = adj[j]; break; }
                            }
                            if (next == UINT32_MAX) break;
                            {
                                uint32_t vb = 0;
                                add_vertex(&s, input->vertices[next], &vb);
                                add_tri(&s, centroid_idx, cur, vb);
                            }
                            prev = cur;
                            cur = next;
                        }
                    }
                    free(used);
                }
            }
            free(count); free(off); free(adj);
        }
        free(edges);
    }

    if (s.icount) {
        float eps = 1e-5f;
        rc = geomlib_soup_to_mesh(s.v, s.count, s.idx, s.icount, eps, output);
    } else {
        geomlib_trimesh_init(output, 0, 0);
        rc = GEOMLIB_MESH_SUCCESS;
    }
    free(s.v); free(s.idx);
    free(centers); free(assign); free(seeds); free(points); free(normals);
    return rc;
}
