/*
 * geomlib - Robust topology repair (see geomlib_repair.h)
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_repair.h"
#include "geomlib/geomlib_topology.h"
#include "geomlib/geomlib_healing.h"

#include <stdlib.h>
#include <string.h>

typedef struct { uint32_t a, b; } edge_key;

static int cmp_edge_key(const void *pa, const void *pb)
{
    const edge_key *x = (const edge_key *)pa;
    const edge_key *y = (const edge_key *)pb;
    if (x->a < y->a) return -1;
    if (x->a > y->a) return 1;
    if (x->b < y->b) return -1;
    if (x->b > y->b) return 1;
    return 0;
}

int geomlib_repair_nonmanifold_edges(geomlib_trimesh *mesh, size_t *detached_out)
{
    size_t detached = 0;
    size_t *counts;
    edge_key *keys;
    size_t t, k, pass;

    if (!mesh || !mesh->indices) return GEOMLIB_MESH_ERROR_INVALID;
    counts = (size_t *)calloc(mesh->triangle_count, sizeof(size_t));
    keys = (edge_key *)malloc(3 * mesh->triangle_count * sizeof(edge_key));
    if (!counts || !keys) {
        free(counts); free(keys);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (pass = 0; pass < mesh->triangle_count; ++pass) {
        size_t ne = 0;
        int found = 0;
        memset(counts, 0, mesh->triangle_count * sizeof(size_t));
        for (t = 0; t < mesh->triangle_count; ++t) {
            uint32_t v[3] = { mesh->indices[t*3], mesh->indices[t*3+1], mesh->indices[t*3+2] };
            for (k = 0; k < 3; ++k) {
                keys[ne].a = v[k] < v[(k+1)%3] ? v[k] : v[(k+1)%3];
                keys[ne].b = v[k] < v[(k+1)%3] ? v[(k+1)%3] : v[k];
                ++ne;
            }
        }
        qsort(keys, ne, sizeof(edge_key), cmp_edge_key);
        for (t = 0; t < ne;) {
            size_t j = t + 1;
            while (j < ne && cmp_edge_key(&keys[t], &keys[j]) == 0) ++j;
            if (j - t > 2) { found = 1; break; }
            t = j;
        }
        if (!found) break;

        /* keys[t] is the over-used edge. Detach the last triangle incident to
         * it by duplicating all three of its vertices. */
        {
            uint32_t a = keys[t].a, b = keys[t].b;
            size_t victim = (size_t)-1;
            for (k = 0; k < mesh->triangle_count; ++k) {
                uint32_t v0 = mesh->indices[k*3], v1 = mesh->indices[k*3+1],
                         v2 = mesh->indices[k*3+2];
                if ((v0 == a && v1 == b) || (v0 == b && v1 == a) ||
                    (v1 == a && v2 == b) || (v1 == b && v2 == a) ||
                    (v2 == a && v0 == b) || (v2 == b && v0 == a)) {
                    victim = k;
                }
            }
            if (victim == (size_t)-1) break;
            {
                uint32_t old[3] = { mesh->indices[victim*3], mesh->indices[victim*3+1],
                                    mesh->indices[victim*3+2] };
                uint32_t fresh[3];
                for (k = 0; k < 3; ++k) {
                    fresh[k] = (uint32_t)mesh->vertex_count;
                    geomlib_trimesh_add_vertex(mesh, mesh->vertices[old[k]], NULL, NULL);
                }
                for (k = 0; k < 3; ++k) mesh->indices[victim*3+k] = fresh[k];
            }
            ++detached;
        }
    }

    free(counts); free(keys);
    /* Detaching triangles can leave vertices with disconnected fans. */
    if (detached > 0) geomlib_split_nonmanifold_vertices(mesh);
    if (detached_out) *detached_out = detached;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_soup_to_mesh(const geomlib_vec3 *points, size_t point_count,
                         const uint32_t *indices, size_t index_count,
                         float weld_eps, geomlib_trimesh *out)
{
    uint32_t *map, *representative;
    geomlib_trimesh mesh;
    size_t i, vcount = 0;
    float eps2 = weld_eps * weld_eps;

    if (!points || !indices || !out || (index_count % 3) != 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    map = (uint32_t *)malloc(point_count * sizeof(uint32_t));
    representative = (uint32_t *)malloc(point_count * sizeof(uint32_t));
    if (!map || !representative) {
        free(map); free(representative);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    /* Weld by comparing each point to the representative point of each
     * unique vertex (O(n^2); suitable for the typical soup sizes here). */
    for (i = 0; i < point_count; ++i) {
        size_t j;
        map[i] = UINT32_MAX;
        for (j = 0; j < vcount; ++j) {
            geomlib_vec3 d = geomlib_vec3_sub(points[i], points[representative[j]]);
            if (geomlib_vec3_dot(d, d) <= eps2) { map[i] = (uint32_t)j; break; }
        }
        if (map[i] == UINT32_MAX) {
            map[i] = (uint32_t)vcount;
            representative[vcount++] = (uint32_t)i;
        }
    }

    if (geomlib_trimesh_init(&mesh, vcount, index_count / 3) != GEOMLIB_MESH_SUCCESS) {
        free(map);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < vcount; ++i)
        geomlib_trimesh_add_vertex(&mesh, points[representative[i]], NULL, NULL);
    {
        /* Skip duplicate faces regardless of winding (sorted vertex triple). */
        edge_key *faces = (edge_key *)malloc((index_count / 3 + 1) * sizeof(edge_key) * 2);
        size_t nfaces = 0;
        if (!faces) {
            geomlib_trimesh_free(&mesh);
            free(map); free(representative);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        for (i = 0; i + 2 < index_count; i += 3) {
            uint32_t a = map[indices[i]], b = map[indices[i+1]], c = map[indices[i+2]];
            uint32_t lo, mid, hi, tmp;
            size_t j;
            int duplicate = 0;
            if (a == b || b == c || a == c) continue;
            lo = a; mid = b; hi = c;
            if (lo > mid) { tmp = lo; lo = mid; mid = tmp; }
            if (mid > hi) { tmp = mid; mid = hi; hi = tmp; }
            if (lo > mid) { tmp = lo; lo = mid; mid = tmp; }
            for (j = 0; j < nfaces; ++j) {
                if (faces[j*2].a == lo && faces[j*2].b == mid &&
                    faces[j*2+1].a == hi) { duplicate = 1; break; }
            }
            if (duplicate) continue;
            faces[nfaces*2].a = lo;
            faces[nfaces*2].b = mid;
            faces[nfaces*2+1].a = hi;
            faces[nfaces*2+1].b = 0;
            ++nfaces;
            geomlib_trimesh_add_triangle(&mesh, a, b, c);
        }
        free(faces);
    }
    free(map);
    free(representative);

    geomlib_remove_duplicate_faces(&mesh);
    geomlib_orient_consistently(&mesh);
    geomlib_remove_isolated_vertices(&mesh);
    *out = mesh;
    return GEOMLIB_MESH_SUCCESS;
}
