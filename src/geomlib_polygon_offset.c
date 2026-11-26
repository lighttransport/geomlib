/*
 * geomlib - 2D polygon offsetting (see geomlib_polygon_offset.h)
 *
 * Clean-room implementation of line-offset polygon offsetting: each edge is
 * translated along its normal by the requested distance, consecutive offset
 * lines are intersected to form the new vertices, and over-long miters are
 * clamped. Self-crossing offsets are split into simple loops.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polygon_offset.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_offset_result geomlib_offset_result_create(void)
{
    geomlib_offset_result r;
    memset(&r, 0, sizeof(r));
    return r;
}

void geomlib_offset_result_free(geomlib_offset_result *r)
{
    size_t i;
    if (!r) return;
    for (i = 0; i < r->count; ++i) geomlib_polygon2d_free(&r->polygons[i]);
    free(r->polygons);
    memset(r, 0, sizeof(*r));
}

typedef struct { geomlib_vec2 p, d; } line;

static geomlib_vec2 offset_vertex(geomlib_vec2 p, geomlib_vec2 n, float dist)
{
    return geomlib_vec2_add(p, geomlib_vec2_mul(n, dist));
}

/* Intersection of two infinite lines; falls back to the midpoint. */
static geomlib_vec2 line_intersection(line a, line b, geomlib_vec2 fallback)
{
    double den = a.d.x*b.d.y - a.d.y*b.d.x;
    if (fabs(den) < 1e-12) return fallback;
    double dx = b.p.x - a.p.x;
    double dy = b.p.y - a.p.y;
    double t = (dx*b.d.y - dy*b.d.x) / den;
    return geomlib_vec2_add(a.p, geomlib_vec2_mul(a.d, (float)t));
}

static int segment_intersection(geomlib_vec2 a, geomlib_vec2 b,
                                geomlib_vec2 c, geomlib_vec2 d,
                                geomlib_vec2 *out)
{
    double r0 = b.x-a.x, r1 = b.y-a.y;
    double s0 = d.x-c.x, s1 = d.y-c.y;
    double den = r0*s1 - r1*s0;
    double t, u;
    if (fabs(den) < 1e-12) return 0;
    t = ((c.x-a.x)*s1 - (c.y-a.y)*s0)/den;
    u = ((c.x-a.x)*r1 - (c.y-a.y)*r0)/den;
    if (t < 1e-9 || t > 1.0-1e-9 || u < 1e-9 || u > 1.0-1e-9) return 0;
    out->x = (float)(a.x + t*r0);
    out->y = (float)(a.y + t*r1);
    return 1;
}

/* Split a possibly self-crossing loop at intersections into simple loops. */
static int split_loops(geomlib_vec2 *v, size_t n, geomlib_offset_result *out)
{
    size_t start = 0;
    int rc = GEOMLIB_MESH_SUCCESS;
    if (n < 3) return GEOMLIB_MESH_SUCCESS;
    while (start + 2 < n) {
        size_t end = start + 1;
        int found = 0;
        for (; end + 1 < n; ++end) {
            size_t i;
            for (i = start; i < end; ++i) {
                geomlib_vec2 hit;
                if (segment_intersection(v[i], v[i+1], v[end], v[end+1], &hit)) {
                    /* Keep the loop up to the crossing and restart there. */
                    size_t len = end - start + 1;
                    if (len >= 3) {
                        geomlib_polygon2d part = {0};
                        size_t k;
                        if (geomlib_polygon2d_init(&part, len) != GEOMLIB_MESH_SUCCESS)
                            return GEOMLIB_MESH_ERROR_MEMORY;
                        for (k = 0; k < len; ++k)
                            geomlib_polygon2d_add_vertex(&part, v[start+k]);
                        {
                            geomlib_polygon2d *np = (geomlib_polygon2d *)realloc(
                                out->polygons, (out->count+1)*sizeof(geomlib_polygon2d));
                            if (!np) { geomlib_polygon2d_free(&part); return GEOMLIB_MESH_ERROR_MEMORY; }
                            out->polygons = np;
                            out->polygons[out->count++] = part;
                        }
                    }
                    /* Replace the remaining edge endpoints with the crossing
                     * point to continue from a clean loop. */
                    v[end] = hit;
                    start = end;
                    found = 1;
                    break;
                }
            }
            if (found) break;
        }
        if (!found) {
            size_t len = n - start;
            if (len >= 3) {
                geomlib_polygon2d part = {0};
                size_t k;
                if (geomlib_polygon2d_init(&part, len) != GEOMLIB_MESH_SUCCESS)
                    return GEOMLIB_MESH_ERROR_MEMORY;
                for (k = 0; k < len; ++k)
                    geomlib_polygon2d_add_vertex(&part, v[start+k]);
                {
                    geomlib_polygon2d *np = (geomlib_polygon2d *)realloc(
                        out->polygons, (out->count+1)*sizeof(geomlib_polygon2d));
                    if (!np) { geomlib_polygon2d_free(&part); return GEOMLIB_MESH_ERROR_MEMORY; }
                    out->polygons = np;
                    out->polygons[out->count++] = part;
                }
            }
            break;
        }
    }
    return rc;
}

int geomlib_offset_polygon(const geomlib_polygon2d *polygon,
                           float distance,
                           int arc_segments,
                           geomlib_offset_result *out)
{
    size_t n, i;
    line *lines;
    geomlib_vec2 *result;
    float signed_dist;
    (void)arc_segments;

    if (!polygon || !out || !polygon->vertices || polygon->vertex_count < 3 ||
        !isfinite(distance))
        return GEOMLIB_MESH_ERROR_INVALID;
    n = polygon->vertex_count;
    lines = (line *)malloc(n*sizeof(line));
    result = (geomlib_vec2 *)malloc(n*sizeof(geomlib_vec2));
    if (!lines || !result) {
        free(lines); free(result);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    /* For a CCW polygon the left normal points inward; positive distance
     * offsets outward, hence the sign flip. */
    signed_dist = -distance;
    for (i = 0; i < n; ++i) {
        geomlib_vec2 p = polygon->vertices[i];
        geomlib_vec2 q = polygon->vertices[(i+1)%n];
        geomlib_vec2 e = geomlib_vec2_sub(q, p);
        float len = sqrtf(e.x*e.x + e.y*e.y);
        geomlib_vec2 inward;
        if (len < 1e-9f) {
            lines[i].p = p;
            lines[i].d = geomlib_vec2_create(1,0);
            continue;
        }
        inward = geomlib_vec2_create(-e.y/len, e.x/len);
        lines[i].p = offset_vertex(p, inward, signed_dist);
        lines[i].d = geomlib_vec2_create(e.x/len, e.y/len);
    }
    for (i = 0; i < n; ++i) {
        geomlib_vec2 fallback = geomlib_vec2_mul(
            geomlib_vec2_add(lines[i].p, lines[(i+1)%n].p), 0.5f);
        result[i] = line_intersection(lines[(i+n-1)%n], lines[i], fallback);
    }

    /* Clamp long miters so spikes cannot grow unbounded. */
    {
        float max_miter = fabsf(distance)*4.0f + 1e-3f;
        for (i = 0; i < n; ++i) {
            geomlib_vec2 d = geomlib_vec2_sub(result[i], polygon->vertices[i]);
            float m = sqrtf(d.x*d.x + d.y*d.y);
            if (m > max_miter) {
                result[i] = geomlib_vec2_add(polygon->vertices[i],
                    geomlib_vec2_mul(d, max_miter/m));
            }
        }
    }

    {
        int rc = split_loops(result, n, out);
        free(lines); free(result);
        return rc;
    }
}
