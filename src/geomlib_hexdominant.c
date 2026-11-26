/*
 * geomlib - Hex-dominant meshing (see geomlib_hexdominant.h)
 *
 * Clean-room implementation of the 5-tet-to-hex merge: find edges shared by
 * exactly five tetrahedra whose union has eight distinct vertices, order the
 * six ring vertices around the diagonal by angle, and emit a hex cell.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_hexdominant.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_hexmesh geomlib_hexmesh_create(void)
{
    geomlib_hexmesh m;
    memset(&m, 0, sizeof(m));
    return m;
}

void geomlib_hexmesh_free(geomlib_hexmesh *m)
{
    if (!m) return;
    free(m->hexes);
    free(m->tets);
    memset(m, 0, sizeof(*m));
}

typedef struct { uint32_t a, b; } ekey;

static int cmp_ekey(const void *pa, const void *pb)
{
    const ekey *x=(const ekey*)pa, *y=(const ekey*)pb;
    if (x->a<y->a) return -1;
    if (x->a>y->a) return 1;
    if (x->b<y->b) return -1;
    if (x->b>y->b) return 1;
    return 0;
}

static int find_vertex(const uint32_t *set, size_t n, uint32_t v)
{
    size_t i;
    for (i = 0; i < n; ++i) if (set[i] == v) return 1;
    return 0;
}

int geomlib_hex_dominant(const geomlib_tetramesh *input,
                         geomlib_hexmesh *out)
{
    ekey *edges;
    size_t *tets_of_edge, ne = 0, i, k, nhex = 0, ntet = 0;
    uint8_t *consumed;
    geomlib_hexmesh result = geomlib_hexmesh_create();

    if (!input || !out || !input->tetrahedra || input->tetrahedron_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    edges = (ekey *)malloc(6*input->tetrahedron_count*sizeof(ekey));
    tets_of_edge = (size_t *)malloc(6*input->tetrahedron_count*sizeof(size_t));
    consumed = (uint8_t *)calloc(input->tetrahedron_count, sizeof(uint8_t));
    if (!edges || !tets_of_edge || !consumed) {
        free(edges); free(tets_of_edge); free(consumed);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < input->tetrahedron_count; ++i) {
        uint32_t v[4];
        for (k = 0; k < 4; ++k) v[k] = input->tetrahedra[i].vertices[k];
        for (k = 0; k < 6; ++k) {
            static const int pairs[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
            uint32_t a=v[pairs[k][0]], b=v[pairs[k][1]];
            edges[ne].a = a<b?a:b;
            edges[ne].b = a<b?b:a;
            tets_of_edge[ne] = i;
            ne++;
        }
    }
    {
        /* Sort edges and their tet ids together (insertion sort). */
        for (i = 1; i < ne; ++i) {
            ekey ke = edges[i];
            size_t kt = tets_of_edge[i];
            size_t j = i;
            while (j > 0 && cmp_ekey(&edges[j-1], &ke) > 0) {
                edges[j] = edges[j-1];
                tets_of_edge[j] = tets_of_edge[j-1];
                --j;
            }
            edges[j] = ke;
            tets_of_edge[j] = kt;
        }
    }

    result.hexes = (uint32_t *)malloc((input->tetrahedron_count/5+1)*8*sizeof(uint32_t));
    if (!result.hexes) { free(edges); free(tets_of_edge); free(consumed); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (i = 0; i < ne;) {
        size_t j = i+1, count, nv = 0, q, r;
        uint32_t verts[16];
        uint32_t ring[6], diagonal[2];
        while (j < ne && cmp_ekey(&edges[i], &edges[j]) == 0) ++j;
        if (j-i != 5) { i = j; continue; }
        diagonal[0] = edges[i].a;
        diagonal[1] = edges[i].b;
        count = 0;
        for (q = i; q < j; ++q) {
            uint32_t t = (uint32_t)tets_of_edge[q];
            if (consumed[t]) { count = 0; break; }
            for (r = 0; r < 4; ++r) {
                uint32_t v = input->tetrahedra[t].vertices[r];
                if (!find_vertex(verts, nv, v)) {
                    if (nv < 16) verts[nv++] = v;
                }
            }
            count++;
        }
        if (count != 5 || nv != 8) { i = j; continue; }
        /* Six ring vertices: all except the two diagonal endpoints. */
        {
            size_t nr = 0;
            for (q = 0; q < nv && nr < 6; ++q)
                if (verts[q] != diagonal[0] && verts[q] != diagonal[1])
                    ring[nr++] = verts[q];
            if (nr != 6) { i = j; continue; }
        }
        /* Sort ring vertices by angle around the diagonal. */
        {
            geomlib_vec3 axis = geomlib_vec3_sub(input->vertices[diagonal[1]],
                                                 input->vertices[diagonal[0]]);
            geomlib_vec3 t1 = geomlib_vec3_sub(input->vertices[ring[0]],
                                               input->vertices[diagonal[0]]);
            geomlib_vec3 t2;
            double center[2] = {0,0};
            geomlib_vec3 c = geomlib_vec3_mul(
                geomlib_vec3_add(input->vertices[diagonal[0]],
                                 input->vertices[diagonal[1]]), 0.5f);
            for (q = 0; q < 6; ++q) {
                geomlib_vec3 d = geomlib_vec3_sub(input->vertices[ring[q]], c);
                center[0] += d.x; center[1] += d.y;
            }
            (void)center;
            axis = geomlib_vec3_normalize(axis);
            t1 = geomlib_vec3_sub(t1, geomlib_vec3_mul(axis, geomlib_vec3_dot(axis,t1)));
            if (geomlib_vec3_length(t1) < 1e-8f) t1 = geomlib_vec3_create(1,0,0);
            t1 = geomlib_vec3_normalize(t1);
            t2 = geomlib_vec3_cross(axis, t1);
            for (q = 1; q < 6; ++q) {
                uint32_t cur = ring[q];
                geomlib_vec3 d = geomlib_vec3_sub(input->vertices[cur], c);
                double ang = atan2(geomlib_vec3_dot(d,t2), geomlib_vec3_dot(d,t1));
                size_t p = q;
                while (p > 0) {
                    uint32_t prev = ring[p-1];
                    geomlib_vec3 dp = geomlib_vec3_sub(input->vertices[prev], c);
                    double ap = atan2(geomlib_vec3_dot(dp,t2), geomlib_vec3_dot(dp,t1));
                    if (ap <= ang) break;
                    ring[p] = prev;
                    --p;
                }
                ring[p] = cur;
            }
        }
        for (q = 0; q < j-i; ++q) consumed[tets_of_edge[i+q]] = 1;
        /* Hex: four ring vertices around one diagonal endpoint, then four
         * around the other (ordering is topological, not geometric). */
        for (q = 0; q < 4; ++q) result.hexes[nhex*8+q] = ring[q];
        result.hexes[nhex*8+4] = diagonal[0];
        for (q = 4; q < 6; ++q) result.hexes[nhex*8+1+q] = ring[q];
        result.hexes[nhex*8+7] = diagonal[1];
        nhex++;
        i = j;
    }
    free(edges);
    free(tets_of_edge);

    result.tets = (uint32_t *)malloc(input->tetrahedron_count*4*sizeof(uint32_t));
    if (!result.tets) { free(consumed); geomlib_hexmesh_free(&result); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (i = 0; i < input->tetrahedron_count; ++i) {
        if (consumed[i]) continue;
        for (k = 0; k < 4; ++k)
            result.tets[ntet*4+k] = input->tetrahedra[i].vertices[k];
        ntet++;
    }
    free(consumed);
    result.hex_count = nhex;
    result.tet_count = ntet;
    result.vertices = input->vertices;
    result.vertex_count = input->vertex_count;
    *out = result;
    return GEOMLIB_MESH_SUCCESS;
}
