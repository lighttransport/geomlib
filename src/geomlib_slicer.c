/*
 * geomlib - Mesh slicing (see geomlib_slicer.h)
 *
 * Clean-room implementation of plane/mesh slicing: each triangle is clipped
 * by the plane, producing 0 or 2 crossing points; segments are chained into
 * polylines by matching endpoints with a welding tolerance.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_slicer.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_slice_result geomlib_slice_result_create(void)
{
    geomlib_slice_result r;
    memset(&r, 0, sizeof(r));
    return r;
}

void geomlib_slice_result_free(geomlib_slice_result *r)
{
    if (!r) return;
    free(r->points);
    free(r->starts);
    free(r->counts);
    memset(r, 0, sizeof(*r));
}

typedef struct { geomlib_vec3 a, b; } segment;

static float plane_value(geomlib_vec3 p, geomlib_vec3 point, geomlib_vec3 normal)
{
    return geomlib_vec3_dot(normal, geomlib_vec3_sub(p, point));
}

static int same_point(geomlib_vec3 a, geomlib_vec3 b, float eps)
{
    geomlib_vec3 d = geomlib_vec3_sub(a, b);
    return geomlib_vec3_dot(d, d) <= eps*eps;
}

int geomlib_slice_mesh(const geomlib_trimesh *mesh,
                       geomlib_vec3 point,
                       geomlib_vec3 normal,
                       float weld_eps,
                       geomlib_slice_result *out)
{
    segment *segments;
    size_t nseg = 0, i, j, k;
    int *used;

    if (!mesh || !out || !mesh->vertices || !mesh->indices ||
        mesh->triangle_count == 0 || geomlib_vec3_length(normal) < 1e-10f)
        return GEOMLIB_MESH_ERROR_INVALID;
    if (!(weld_eps > 0.0f)) weld_eps = 1e-4f;

    segments = (segment *)malloc(mesh->triangle_count*sizeof(segment));
    used = (int *)calloc(mesh->triangle_count, sizeof(int));
    if (!segments || !used) {
        free(segments); free(used);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (i = 0; i < mesh->triangle_count; ++i) {
        uint32_t v[3] = { mesh->indices[i*3], mesh->indices[i*3+1], mesh->indices[i*3+2] };
        float d[3];
        geomlib_vec3 crossings[4];
        size_t nc = 0;
        for (k = 0; k < 3; ++k) d[k] = plane_value(mesh->vertices[v[k]], point, normal);
        for (k = 0; k < 3; ++k) {
            size_t next = (k+1)%3;
            if ((d[k] < 0.0f) == (d[next] < 0.0f)) continue;
            if (fabsf(d[k]-d[next]) < 1e-12f) continue;
            float t = d[k] / (d[k]-d[next]);
            crossings[nc++] = geomlib_vec3_add(
                geomlib_vec3_mul(mesh->vertices[v[k]], 1.0f-t),
                geomlib_vec3_mul(mesh->vertices[v[next]], t));
        }
        if (nc == 2) { segments[nseg].a = crossings[0]; segments[nseg].b = crossings[1]; nseg++; }
    }
    free(used);
    if (!nseg) {
        free(segments);
        memset(out, 0, sizeof(*out));
        return GEOMLIB_MESH_SUCCESS;
    }

    /* Chain segments into polylines. */
    {
        int *consumed = (int *)calloc(nseg, sizeof(int));
        geomlib_vec3 *points = (geomlib_vec3 *)malloc(2*nseg*sizeof(geomlib_vec3));
        uint32_t *starts = (uint32_t *)malloc(nseg*sizeof(uint32_t));
        uint32_t *counts = (uint32_t *)calloc(nseg, sizeof(uint32_t));
        size_t np = 0, nlines = 0;
        if (!consumed || !points || !starts || !counts) {
            free(consumed); free(points); free(starts); free(counts); free(segments);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        for (i = 0; i < nseg; ++i) {
            geomlib_vec3 cur_a, cur_b;
            if (consumed[i]) continue;
            consumed[i] = 1;
            cur_a = segments[i].a;
            cur_b = segments[i].b;
            starts[nlines] = (uint32_t)np;
            points[np++] = cur_a;
            points[np++] = cur_b;
            counts[nlines] = 2;
            /* Extend at the tail until no matching segment remains. */
            for (;;) {
                int found = -1;
                int at_b = 0;
                for (j = 0; j < nseg; ++j) {
                    if (consumed[j]) continue;
                    if (same_point(segments[j].a, cur_b, weld_eps)) { found = (int)j; at_b = 0; break; }
                    if (same_point(segments[j].b, cur_b, weld_eps)) { found = (int)j; at_b = 1; break; }
                }
                if (found < 0) break;
                consumed[found] = 1;
                cur_b = at_b ? segments[found].a : segments[found].b;
                if (!same_point(cur_b, points[np-1], weld_eps)) {
                    points[np++] = cur_b;
                    counts[nlines]++;
                }
            }
            /* Prepend at the head. */
            for (;;) {
                int found = -1;
                int at_b = 0;
                for (j = 0; j < nseg; ++j) {
                    if (consumed[j]) continue;
                    if (same_point(segments[j].b, cur_a, weld_eps)) { found = (int)j; at_b = 0; break; }
                    if (same_point(segments[j].a, cur_a, weld_eps)) { found = (int)j; at_b = 1; break; }
                }
                if (found < 0) break;
                consumed[found] = 1;
                cur_a = at_b ? segments[found].b : segments[found].a;
                {
                    size_t pos = starts[nlines];
                    size_t cnt = counts[nlines];
                    for (j = 0; j < cnt; ++j) points[pos+cnt-j] = points[pos+cnt-j-1];
                    points[pos] = cur_a;
                    counts[nlines]++;
                }
            }
            nlines++;
        }
        out->points = points;
        out->point_count = np;
        out->starts = starts;
        out->counts = counts;
        out->polyline_count = nlines;
        free(consumed);
    }
    free(segments);
    return GEOMLIB_MESH_SUCCESS;
}
