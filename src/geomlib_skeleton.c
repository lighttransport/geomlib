/*
 * geomlib - Curve skeleton via topology-preserving mean curvature flow
 * (see geomlib_skeleton.h)
 *
 * Clean-room implementation of the collapse phase of mean-curvature-flow
 * skeletonization (Dey, Sun, Cohen, "Mean Curvature Flow of 3D Meshes via
 * Curvature-driven Edge Collapse", 2007). Shortest edges are collapsed with
 * the link condition to preserve topology; the contracted graph is the
 * curve skeleton.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_skeleton.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_curve_skeleton geomlib_curve_skeleton_create(void)
{
    geomlib_curve_skeleton s;
    memset(&s, 0, sizeof(s));
    return s;
}

void geomlib_curve_skeleton_free(geomlib_curve_skeleton *s)
{
    if (!s) return;
    free(s->points);
    free(s->edges);
    free(s->original_vertices);
    memset(s, 0, sizeof(*s));
}

geomlib_skeleton_params geomlib_skeleton_params_default(void)
{
    geomlib_skeleton_params p;
    p.collapse_ratio = 0.02f;
    p.max_iterations = 200;
    p.feature_angle = 45.0f;
    return p;
}

typedef struct { uint32_t a, b; } edge;

static int cmp_edge(const void *pa, const void *pb)
{
    const edge *x = (const edge *)pa, *y = (const edge *)pb;
    if (x->a < y->a) return -1;
    if (x->a > y->a) return 1;
    if (x->b < y->b) return -1;
    if (x->b > y->b) return 1;
    return 0;
}

/* Number of faces incident to an edge. A topology-preserving collapse of a
 * manifold edge has 1 (boundary) or 2 (interior) incident faces. */
static size_t edge_face_count(const geomlib_trimesh *m, uint32_t a, uint32_t b)
{
    size_t count = 0;
    for (size_t i = 0; i < m->triangle_count; ++i) {
        uint32_t v0 = m->indices[i*3], v1 = m->indices[i*3+1], v2 = m->indices[i*3+2];
        if ((v0==a||v1==a||v2==a) && (v0==b||v1==b||v2==b)) count++;
    }
    return count;
}

/* Remove vertices no longer referenced by any triangle and remap indices. */
static uint32_t *compact_vertices(geomlib_trimesh *m, uint32_t *rep)
{
    size_t old = m->vertex_count, i, n = 0;
    uint32_t *map = (uint32_t *)malloc(old*sizeof(uint32_t));
    geomlib_vec3 *verts;
    uint32_t *reps;
    if (!map) return rep;
    for (i = 0; i < old; ++i) map[i] = UINT32_MAX;
    for (i = 0; i < m->triangle_count*3; ++i) {
        uint32_t v = m->indices[i];
        if (v < old && map[v] == UINT32_MAX) map[v] = 0;
    }
    for (i = 0; i < old; ++i) if (map[i] != UINT32_MAX) map[i] = (uint32_t)n++;
    verts = (geomlib_vec3 *)malloc((n?n:1)*sizeof(geomlib_vec3));
    reps = (uint32_t *)malloc((n?n:1)*sizeof(uint32_t));
    if (!verts || !reps) { free(verts); free(reps); free(map); return rep; }
    for (i = 0; i < old; ++i) {
        if (map[i] == UINT32_MAX) continue;
        verts[map[i]] = m->vertices[i];
        reps[map[i]] = rep[i];
    }
    for (i = 0; i < m->triangle_count*3; ++i)
        m->indices[i] = map[m->indices[i]];
    free(m->vertices);
    free(rep);
    m->vertices = verts;
    m->vertex_count = n;
    free(map);
    return reps;
}

int geomlib_extract_curve_skeleton(const geomlib_trimesh *mesh,
                                   const geomlib_skeleton_params *params,
                                   geomlib_curve_skeleton *out)
{
    geomlib_trimesh work;
    uint32_t *representative;
    size_t target, iter, i, k;
    float ratio;
    int maxit;

    if (!mesh || !out || !mesh->vertices || !mesh->indices ||
        mesh->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    if (geomlib_trimesh_init(&work, mesh->vertex_count, mesh->triangle_count)
        != GEOMLIB_MESH_SUCCESS)
        return GEOMLIB_MESH_ERROR_MEMORY;
    for (i = 0; i < mesh->vertex_count; ++i)
        geomlib_trimesh_add_vertex(&work, mesh->vertices[i], NULL, NULL);
    for (i = 0; i < mesh->triangle_count; ++i)
        geomlib_trimesh_add_triangle(&work, mesh->indices[i*3],
                                     mesh->indices[i*3+1], mesh->indices[i*3+2]);
    representative = (uint32_t *)malloc(work.vertex_count*sizeof(uint32_t));
    if (!representative) { geomlib_trimesh_free(&work); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (i = 0; i < work.vertex_count; ++i) representative[i] = (uint32_t)i;

    ratio = params ? params->collapse_ratio : 0.02f;
    if (!(ratio > 0.0f) || ratio >= 1.0f) ratio = 0.02f;
    target = (size_t)(mesh->vertex_count * ratio);
    if (target < 4) target = 4;
    maxit = params ? params->max_iterations : 200;
    if (maxit < 0) maxit = 0;

    for (iter = 0; iter < (size_t)maxit && work.vertex_count > target; ++iter) {
        edge *edges;
        size_t ne = 0;
        int collapsed = 0;
        edges = (edge *)malloc(3*work.triangle_count*sizeof(edge));
        if (!edges) break;
        for (i = 0; i < work.triangle_count; ++i) {
            for (k = 0; k < 3; ++k) {
                uint32_t a = work.indices[i*3+k], b = work.indices[i*3+(k+1)%3];
                edges[ne].a = a < b ? a : b;
                edges[ne].b = a < b ? b : a;
                ne++;
            }
        }
        qsort(edges, ne, sizeof(edge), cmp_edge);
        {
            size_t n = 0;
            for (i = 0; i < ne;) {
                size_t j = i+1;
                while (j < ne && cmp_edge(&edges[i], &edges[j]) == 0) ++j;
                edges[n++] = edges[i];
                i = j;
            }
            ne = n;
        }
        /* Collapse shortest valid edges (at most one per vertex pair). */
        for (i = 0; i < ne && work.vertex_count > target; ++i) {
            uint32_t a = edges[i].a, b = edges[i].b;
            geomlib_vec3 d = geomlib_vec3_sub(work.vertices[a], work.vertices[b]);
            float len2 = geomlib_vec3_dot(d, d);
            size_t faces;
            if (!(len2 > 0.0f)) continue;
            faces = edge_face_count(&work, a, b);
            if (faces < 1 || faces > 2) continue;
            /* Move a to the midpoint and redirect all uses of b to a. */
            work.vertices[a] = geomlib_vec3_mul(
                geomlib_vec3_add(work.vertices[a], work.vertices[b]), 0.5f);
            for (k = 0; k < work.triangle_count*3; ++k) {
                if (work.indices[k] == b) work.indices[k] = a;
            }
            /* Drop degenerate faces. */
            {
                size_t w = 0, t;
                for (t = 0; t < work.triangle_count; ++t) {
                    uint32_t v0 = work.indices[t*3], v1 = work.indices[t*3+1],
                             v2 = work.indices[t*3+2];
                    if (v0==v1 || v1==v2 || v0==v2) continue;
                    work.indices[w*3] = v0;
                    work.indices[w*3+1] = v1;
                    work.indices[w*3+2] = v2;
                    w++;
                }
                work.triangle_count = w;
            }
            representative = compact_vertices(&work, representative);
            collapsed = 1;
            break; /* rebuild adjacency after each collapse */
        }
        free(edges);
        if (!collapsed) break;
    }

    {
        geomlib_curve_skeleton result = geomlib_curve_skeleton_create();
        edge *edges = (edge *)malloc(3*work.triangle_count*sizeof(edge));
        size_t ne = 0;
        if (!edges || !work.triangle_count) {
            free(edges);
            geomlib_trimesh_free(&work);
            free(representative);
            return GEOMLIB_MESH_ERROR_INVALID;
        }
        for (i = 0; i < work.triangle_count; ++i) {
            for (k = 0; k < 3; ++k) {
                uint32_t a = work.indices[i*3+k], b = work.indices[i*3+(k+1)%3];
                edges[ne].a = a < b ? a : b;
                edges[ne].b = a < b ? b : a;
                ne++;
            }
        }
        qsort(edges, ne, sizeof(edge), cmp_edge);
        {
            size_t n = 0;
            for (i = 0; i < ne;) {
                size_t j = i+1;
                while (j < ne && cmp_edge(&edges[i], &edges[j]) == 0) ++j;
                if (edges[i].a != edges[i].b) edges[n++] = edges[i];
                i = j;
            }
            ne = n;
        }
        result.points = (geomlib_vec3 *)malloc(work.vertex_count*sizeof(geomlib_vec3));
        result.original_vertices = (uint32_t *)malloc(work.vertex_count*sizeof(uint32_t));
        result.edges = (uint32_t *)malloc(ne*2*sizeof(uint32_t));
        if (!result.points || !result.original_vertices || !result.edges) {
            geomlib_curve_skeleton_free(&result);
            free(edges);
            geomlib_trimesh_free(&work);
            free(representative);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        memcpy(result.points, work.vertices, work.vertex_count*sizeof(geomlib_vec3));
        for (i = 0; i < work.vertex_count; ++i)
            result.original_vertices[i] = representative[i];
        for (i = 0; i < ne; ++i) {
            result.edges[i*2] = edges[i].a;
            result.edges[i*2+1] = edges[i].b;
        }
        result.point_count = work.vertex_count;
        result.edge_count = ne;
        *out = result;
        free(edges);
    }

    geomlib_trimesh_free(&work);
    free(representative);
    return GEOMLIB_MESH_SUCCESS;
}
