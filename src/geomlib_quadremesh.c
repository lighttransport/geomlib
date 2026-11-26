/*
 * geomlib - Quad remeshing (see geomlib_quadremesh.h)
 *
 * Clean-room implementation of greedy triangle pairing for quad remeshing
 * (see e.g. Botsch et al., Polygon Mesh Processing, 2010, chapter on remeshing).
 * Adjacent triangles are merged into a quad when the dihedral angle and the
 * resulting quad angles satisfy quality bounds. An optional iterative
 * planarity relaxation adjusts shared vertices.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_quadremesh.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

geomlib_quadremesh_params geomlib_quadremesh_params_default(void)
{
    geomlib_quadremesh_params p;
    p.max_dihedral = 30.0f;
    p.min_quad_angle = 30.0f;
    p.max_quad_angle = 150.0f;
    p.smoothing_iterations = 0;
    return p;
}

void geomlib_quadmesh_free(geomlib_quadmesh *m)
{
    if (!m) return;
    free(m->vertices);
    free(m->quads);
    free(m->triangles);
    memset(m, 0, sizeof(*m));
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

static float angle_between(geomlib_vec3 a, geomlib_vec3 b, geomlib_vec3 c)
{
    geomlib_vec3 v1 = geomlib_vec3_sub(a, b);
    geomlib_vec3 v2 = geomlib_vec3_sub(c, b);
    float l1 = geomlib_vec3_length(v1), l2 = geomlib_vec3_length(v2);
    if (!(l1 > 1e-10f) || !(l2 > 1e-10f)) return 0.0f;
    float d = geomlib_vec3_dot(v1, v2) / (l1*l2);
    if (d < -1.0f) d = -1.0f;
    if (d > 1.0f) d = 1.0f;
    return acosf(d) * 180.0f / 3.14159265f;
}

int geomlib_remesh_to_quads(const geomlib_trimesh *input,
                            const geomlib_quadremesh_params *params,
                            geomlib_quadmesh *output)
{
    size_t T, E = 0, i, k;
    ekey *edges;
    uint32_t *left, *right;
    size_t *map;
    bool *merged, *used;
    geomlib_vec3 *verts;
    uint32_t *quads, *tris;
    size_t nquads = 0, ntris = 0;
    float cos_limit;

    if (!input || !output || !input->vertices || !input->indices ||
        input->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    T = input->triangle_count;
    edges = (ekey *)malloc(3*T*sizeof(ekey));
    left = (uint32_t *)calloc(3*T, sizeof(uint32_t));
    right = (uint32_t *)calloc(3*T, sizeof(uint32_t));
    map = (size_t *)malloc(3*T*sizeof(size_t));
    if (!edges || !left || !right || !map) {
        free(edges); free(left); free(right); free(map);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < T; ++i) {
        for (k = 0; k < 3; ++k) {
            uint32_t a = input->indices[i*3+k], b = input->indices[i*3+(k+1)%3];
            edges[E].a = a < b ? a : b;
            edges[E].b = a < b ? b : a;
            map[E] = i*3+k;
            E++;
        }
    }
    {
        ekey *sorted = (ekey *)malloc(E*sizeof(ekey));
        size_t *smap = (size_t *)malloc(E*sizeof(size_t));
        size_t *orig = (size_t *)malloc(E*sizeof(size_t));
        if (!sorted || !smap || !orig) {
            free(sorted); free(smap); free(orig);
            free(edges); free(left); free(right); free(map);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        for (i = 0; i < E; ++i) {
            sorted[i] = edges[i];
            smap[i] = map[i];
            orig[i] = i;
        }
        /* Sort keys and their triangle-corner map together (simple insertion
         * sort; E is modest for interactive meshes). */
        for (i = 1; i < E; ++i) {
            ekey ke = sorted[i];
            size_t km = smap[i], ko = orig[i];
            size_t j = i;
            while (j > 0 && cmp_ekey(&sorted[j-1], &ke) > 0) {
                sorted[j] = sorted[j-1];
                smap[j] = smap[j-1];
                orig[j] = orig[j-1];
                --j;
            }
            sorted[j] = ke;
            smap[j] = km;
            orig[j] = ko;
        }
        for (i = 0; i < E;) {
            size_t j = i+1;
            while (j < E && cmp_ekey(&sorted[i], &sorted[j]) == 0) ++j;
            if (j - i == 2) {
                left[orig[i]] = (uint32_t)(smap[i]/3);
                right[orig[i]] = (uint32_t)(smap[i+1]/3);
            }
            i = j;
        }
        free(sorted); free(smap); free(orig);
    }

    merged = (bool *)calloc(T, sizeof(bool));
    used = (bool *)calloc(T, sizeof(bool));
    verts = (geomlib_vec3 *)malloc(input->vertex_count*sizeof(geomlib_vec3));
    quads = (uint32_t *)malloc((T/2+1)*4*sizeof(uint32_t));
    tris = (uint32_t *)malloc(T*3*sizeof(uint32_t));
    if (!merged || !used || !verts || !quads || !tris) {
        free(merged); free(used); free(verts); free(quads); free(tris);
        free(edges); free(left); free(right); free(map);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    memcpy(verts, input->vertices, input->vertex_count*sizeof(geomlib_vec3));
    cos_limit = cosf((params ? params->max_dihedral : 30.0f) * 3.14159265f/180.0f);
    float min_a = params ? params->min_quad_angle : 30.0f;
    float max_a = params ? params->max_quad_angle : 150.0f;

    for (i = 0; i < E; ++i) {
        uint32_t t0 = left[i], t1 = right[i];
        uint32_t q[4];
        geomlib_vec3 n0, n1;
        float min_angle = 180.0f, max_angle = 0.0f;
        if (t0 == UINT32_MAX || t1 == UINT32_MAX || used[t0] || used[t1]) continue;
        if (t0 == t1) continue;
        /* The shared edge endpoints must not be the diagonal of the quad:
         * collect the two opposite vertices of the triangle pair. */
        {
            uint32_t a = edges[i].a, b = edges[i].b;
            uint32_t o0 = UINT32_MAX, o1 = UINT32_MAX;
            for (k = 0; k < 3; ++k) {
                uint32_t v = input->indices[t0*3+k];
                if (v != a && v != b) o0 = v;
                v = input->indices[t1*3+k];
                if (v != a && v != b) o1 = v;
            }
            if (o0 == UINT32_MAX || o1 == UINT32_MAX || o0 == o1) continue;
            q[0] = o0; q[1] = a; q[2] = o1; q[3] = b;
        }
        n0 = geomlib_trimesh_triangle_normal(input, t0);
        n1 = geomlib_trimesh_triangle_normal(input, t1);
        if (geomlib_vec3_dot(n0, n1) < cos_limit) continue;
        for (k = 0; k < 4; ++k) {
            float ang = angle_between(verts[q[(k+3)%4]], verts[q[k]], verts[q[(k+1)%4]]);
            if (ang < min_angle) min_angle = ang;
            if (ang > max_angle) max_angle = ang;
        }
        if (min_angle < min_a || max_angle > max_a) continue;
        for (k = 0; k < 4; ++k) quads[nquads*4+k] = q[k];
        nquads++;
        used[t0] = used[t1] = true;
    }
    for (i = 0; i < T; ++i) {
        if (used[i]) continue;
        for (k = 0; k < 3; ++k) tris[ntris*3+k] = input->indices[i*3+k];
        ntris++;
    }

    /* Optional planarity relaxation of quad vertices (equal weights toward
     * each incident quad's average plane point). */
    {
        int it;
        int relax = params ? params->smoothing_iterations : 0;
        size_t *count = (size_t *)calloc(input->vertex_count, sizeof(size_t));
        if (relax < 0) relax = 0;
        if (relax > 20) relax = 20;
        for (it = 0; it < relax && count; ++it) {
            geomlib_vec3 *acc = (geomlib_vec3 *)calloc(input->vertex_count,
                                                       sizeof(geomlib_vec3));
            if (!acc) break;
            for (i = 0; i < nquads; ++i) {
                geomlib_vec3 c = {0,0,0};
                for (k = 0; k < 4; ++k) {
                    c.x += verts[quads[i*4+k]].x;
                    c.y += verts[quads[i*4+k]].y;
                    c.z += verts[quads[i*4+k]].z;
                }
                c = geomlib_vec3_mul(c, 0.25f);
                for (k = 0; k < 4; ++k) {
                    uint32_t v = quads[i*4+k];
                    acc[v] = geomlib_vec3_add(acc[v], c);
                    count[v]++;
                }
            }
            for (i = 0; i < input->vertex_count; ++i) {
                if (count[i] == 0) continue;
                geomlib_vec3 target = geomlib_vec3_mul(acc[i], 1.0f/(float)count[i]);
                verts[i] = geomlib_vec3_mul(
                    geomlib_vec3_add(verts[i], target), 0.5f);
            }
            free(acc);
        }
        free(count);
    }

    output->vertices = verts;
    output->vertex_count = input->vertex_count;
    output->quads = quads;
    output->quad_count = nquads;
    output->triangles = tris;
    output->triangle_count = ntris;

    free(edges); free(left); free(right); free(map);
    free(merged); free(used);
    return GEOMLIB_MESH_SUCCESS;
}
