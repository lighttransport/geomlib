/*
 * geomlib - Catmull-Clark subdivision (one level)
 *
 * Clean-room implementation of the standard Catmull-Clark rules generalized
 * to polygonal faces; triangular input faces are split into n quads, and each
 * quad is triangulated into two triangles. Boundary vertices use the crease
 * rule. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_subdivision.h"
#include <stdio.h>

#include <stdlib.h>
#include <string.h>

typedef struct {
    uint32_t a, b;
    uint32_t f0, f1;
} cc_edge;

typedef struct {
    uint32_t a, b, face;
} edge_ref;

static int cmp_edge_ref(const void *pa, const void *pb)
{
    const edge_ref *x = (const edge_ref *)pa;
    const edge_ref *y = (const edge_ref *)pb;
    if (x->a < y->a) return -1;
    if (x->a > y->a) return 1;
    if (x->b < y->b) return -1;
    if (x->b > y->b) return 1;
    return 0;
}

static uint32_t find_edge(const cc_edge *edges, size_t n, uint32_t a, uint32_t b)
{
    size_t lo = 0, hi = n;
    if (a > b) { uint32_t t = a; a = b; b = t; }
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (edges[mid].a < a || (edges[mid].a == a && edges[mid].b < b)) lo = mid + 1;
        else hi = mid;
    }
    return (uint32_t)lo;
}

static geomlib_vec3 scale_add(geomlib_vec3 a, geomlib_vec3 b, float s)
{
    return geomlib_vec3_add(a, geomlib_vec3_mul(b, s));
}

int geomlib_subdivide_catmull_clark(const geomlib_trimesh *input,
                                    geomlib_trimesh *output)
{
    size_t V, T, E = 0, i, k;
    edge_ref *refs;
    cc_edge *edges;
    geomlib_vec3 *face_points, *edge_points, *new_vertices;
    size_t *vedge_count, *vface_count, *vface_fill;
    uint32_t *vedges, *vfaces;
    size_t *vedge_off, *vface_off;

    if (!input || !output || !input->vertices || !input->indices ||
        input->vertex_count == 0 || input->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    V = input->vertex_count;
    T = input->triangle_count;

    refs = (edge_ref *)malloc(3 * T * sizeof(edge_ref));
    if (!refs) return GEOMLIB_MESH_ERROR_MEMORY;
    for (i = 0; i < T; ++i) {
        uint32_t v[3] = { input->indices[i*3], input->indices[i*3+1], input->indices[i*3+2] };
        for (k = 0; k < 3; ++k) {
            uint32_t a = v[k], b = v[(k+1)%3];
            refs[E].a = a < b ? a : b;
            refs[E].b = a < b ? b : a;
            refs[E].face = (uint32_t)i;
            ++E;
        }
    }
    qsort(refs, E, sizeof(edge_ref), cmp_edge_ref);
    edges = (cc_edge *)malloc(E * sizeof(cc_edge));
    if (!edges) { free(refs); return GEOMLIB_MESH_ERROR_MEMORY; }
    {
        size_t n = 0;
        for (i = 0; i < E;) {
            size_t j = i + 1;
            while (j < E && refs[i].a == refs[j].a && refs[i].b == refs[j].b) ++j;
            edges[n].a = refs[i].a;
            edges[n].b = refs[i].b;
            edges[n].f0 = refs[i].face;
            /* With two incident faces the second face is refs[i+1].face;
             * more than two is non-manifold and handled as boundary. */
            edges[n].f1 = (j - i == 2) ? refs[i+1].face : UINT32_MAX;
            ++n;
            i = j;
        }
        E = n;
    }
    free(refs);

    face_points = (geomlib_vec3 *)malloc(T * sizeof(geomlib_vec3));
    edge_points = (geomlib_vec3 *)malloc(E * sizeof(geomlib_vec3));
    new_vertices = (geomlib_vec3 *)malloc(V * sizeof(geomlib_vec3));
    vedge_count = (size_t *)calloc(V, sizeof(size_t));
    vface_count = (size_t *)calloc(V, sizeof(size_t));
    if (!face_points || !edge_points || !new_vertices || !vedge_count || !vface_count) {
        free(face_points); free(edge_points); free(new_vertices);
        free(vedge_count); free(vface_count); free(edges);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (i = 0; i < T; ++i) {
        uint32_t a = input->indices[i*3], b = input->indices[i*3+1], c = input->indices[i*3+2];
        face_points[i] = geomlib_vec3_mul(
            geomlib_vec3_add(geomlib_vec3_add(input->vertices[a], input->vertices[b]),
                             input->vertices[c]), 1.0f/3.0f);
        vface_count[a]++; vface_count[b]++; vface_count[c]++;
    }
    for (i = 0; i < E; ++i) {
        vedge_count[edges[i].a]++;
        vedge_count[edges[i].b]++;
    }

    vedge_off = (size_t *)malloc((V+1) * sizeof(size_t));
    vface_off = (size_t *)malloc((V+1) * sizeof(size_t));
    vedges = (uint32_t *)malloc(2 * E * sizeof(uint32_t));
    vfaces = (uint32_t *)malloc(3 * T * sizeof(uint32_t));
    vface_fill = (size_t *)malloc(V * sizeof(size_t));
    if (!vedge_off || !vface_off || !vedges || !vfaces || !vface_fill) {
        free(face_points); free(edge_points); free(new_vertices);
        free(vedge_count); free(vface_count); free(vedge_off); free(vface_off);
        free(vedges); free(vfaces); free(vface_fill); free(edges);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    vedge_off[0] = vface_off[0] = 0;
    for (i = 0; i < V; ++i) {
        vedge_off[i+1] = vedge_off[i] + vedge_count[i];
        vface_off[i+1] = vface_off[i] + vface_count[i];
        vface_fill[i] = 0;
    }
    {
        size_t *cursor = (size_t *)calloc(V, sizeof(size_t));
        if (!cursor) {
            free(face_points); free(edge_points); free(new_vertices);
            free(vedge_count); free(vface_count); free(vedge_off); free(vface_off);
            free(vedges); free(vfaces); free(vface_fill); free(edges);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        for (i = 0; i < E; ++i) {
            vedges[vedge_off[edges[i].a] + cursor[edges[i].a]++] = (uint32_t)i;
            vedges[vedge_off[edges[i].b] + cursor[edges[i].b]++] = (uint32_t)i;
        }
        free(cursor);
    }
    for (i = 0; i < T; ++i) {
        uint32_t v[3] = { input->indices[i*3], input->indices[i*3+1], input->indices[i*3+2] };
        for (k = 0; k < 3; ++k) {
            uint32_t f = vfaces[vface_off[v[k]] + vface_fill[v[k]]++] = (uint32_t)i;
            (void)f;
        }
    }

    /* Edge points. */
    for (i = 0; i < E; ++i) {
        cc_edge *e = &edges[i];
        geomlib_vec3 p = geomlib_vec3_mul(
            geomlib_vec3_add(input->vertices[e->a], input->vertices[e->b]), 0.5f);
        if (e->f1 != UINT32_MAX) {
            p = scale_add(p, face_points[e->f0], 0.25f);
            p = scale_add(p, face_points[e->f1], 0.25f);
        }
        edge_points[i] = p;
    }

    /* Vertex points. */
    for (i = 0; i < V; ++i) {
        size_t ne = vedge_count[i], nf = vface_off[i+1] - vface_off[i];
        geomlib_vec3 q = {0,0,0}, r = {0,0,0};
        size_t j;
        int boundary = 0;
        uint32_t b0 = UINT32_MAX, b1 = UINT32_MAX;
        for (j = 0; j < ne; ++j) {
            uint32_t ei = vedges[vedge_off[i] + j];
            if (ei >= E) return GEOMLIB_MESH_ERROR_INVALID;
            cc_edge *e = &edges[ei];
            uint32_t other = e->a == i ? e->b : e->a;
            geomlib_vec3 mid = geomlib_vec3_mul(
                geomlib_vec3_add(input->vertices[i], input->vertices[other]), 0.5f);
            r = geomlib_vec3_add(r, mid);
            if (e->f1 == UINT32_MAX) {
                boundary = 1;
                if (b0 == UINT32_MAX) b0 = other; else if (b1 == UINT32_MAX) b1 = other;
            }
        }
        for (j = 0; j < nf; ++j) {
            q = geomlib_vec3_add(q, face_points[vfaces[vface_off[i] + j]]);
        }
        if (boundary) {
            geomlib_vec3 v = input->vertices[i];
            if (b0 != UINT32_MAX && b1 != UINT32_MAX)
                new_vertices[i] = geomlib_vec3_mul(
                    geomlib_vec3_add(geomlib_vec3_add(input->vertices[b0],
                                                      geomlib_vec3_mul(v, 6.0f)),
                                     input->vertices[b1]), 1.0f/8.0f);
            else
                new_vertices[i] = v;
        } else if (ne > 0 && nf > 0) {
            float inv = 1.0f / (float)ne;
            geomlib_vec3 vr = geomlib_vec3_mul(input->vertices[i], (float)(ne - 3));
            new_vertices[i] = geomlib_vec3_mul(
                geomlib_vec3_add(geomlib_vec3_add(geomlib_vec3_mul(q, inv),
                                                  geomlib_vec3_mul(r, 2.0f*inv)), vr), inv);
        } else {
            new_vertices[i] = input->vertices[i];
        }
    }

    /* Assemble output vertices and faces. */
    {
        geomlib_trimesh tmp;
        int rc = geomlib_trimesh_init(&tmp, V + E + T, T * 6);
        if (rc != GEOMLIB_MESH_SUCCESS) {
            free(face_points); free(edge_points); free(new_vertices);
            free(vedge_count); free(vface_count); free(vedge_off); free(vface_off);
            free(vedges); free(vfaces); free(vface_fill); free(edges);
            return rc;
        }
        for (i = 0; i < V; ++i) geomlib_trimesh_add_vertex(&tmp, new_vertices[i], NULL, NULL);
        for (i = 0; i < E; ++i) geomlib_trimesh_add_vertex(&tmp, edge_points[i], NULL, NULL);
        for (i = 0; i < T; ++i) geomlib_trimesh_add_vertex(&tmp, face_points[i], NULL, NULL);
        for (i = 0; i < T; ++i) {
            uint32_t corners[3] = { input->indices[i*3], input->indices[i*3+1], input->indices[i*3+2] };
            uint32_t fidx = (uint32_t)(V + E + i);
            for (k = 0; k < 3; ++k) {
                uint32_t cur = corners[k];
                uint32_t next = corners[(k+1)%3];
                uint32_t prev = corners[(k+2)%3];
                uint32_t enext = (uint32_t)(V + find_edge(edges, E, cur, next));
                uint32_t eprev = (uint32_t)(V + find_edge(edges, E, prev, cur));
                geomlib_trimesh_add_triangle(&tmp, cur, enext, fidx);
                geomlib_trimesh_add_triangle(&tmp, cur, fidx, eprev);
            }
        }
        *output = tmp;
    }

    free(face_points); free(edge_points); free(new_vertices);
    free(vedge_count); free(vface_count); free(vedge_off); free(vface_off);
    free(vedges); free(vfaces); free(vface_fill); free(edges);
    return GEOMLIB_MESH_SUCCESS;
}
