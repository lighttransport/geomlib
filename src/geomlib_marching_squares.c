/*
 * geomlib - Marching squares (see geomlib_marching_squares.h)
 *
 * Clean-room implementation of the 16-case marching squares table with a
 * cell-centre saddle decider, followed by endpoint chaining into polylines.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_marching_squares.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_contour_result geomlib_contour_result_create(void)
{
    geomlib_contour_result r;
    memset(&r, 0, sizeof(r));
    return r;
}

void geomlib_contour_result_free(geomlib_contour_result *r)
{
    if (!r) return;
    free(r->points);
    free(r->starts);
    free(r->counts);
    memset(r, 0, sizeof(*r));
}

typedef struct { geomlib_vec2 a, b; } seg;

static int same(geomlib_vec2 a, geomlib_vec2 b, float eps)
{
    geomlib_vec2 d = geomlib_vec2_sub(a,b);
    return geomlib_vec2_dot(d,d) <= eps*eps;
}

static geomlib_vec2 lerp2(geomlib_vec2 a, geomlib_vec2 b, double va, double vb)
{
    double t = (0.0 - va) / (vb - va);
    if (!(t > 0.0 && t < 1.0)) t = 0.5;
    return geomlib_vec2_add(geomlib_vec2_mul(a, (float)(1.0-t)),
                            geomlib_vec2_mul(b, (float)t));
}

int geomlib_marching_squares(const float *values, size_t nx, size_t ny,
                             float h, geomlib_vec2 origin,
                             geomlib_contour_result *out)
{
    seg *segments;
    size_t nseg = 0, x, y, i;
    float eps;

    if (!values || !out || nx < 2 || ny < 2 || !(h > 0.0f))
        return GEOMLIB_MESH_ERROR_INVALID;
    eps = fmaxf(h*1e-4f, 1e-6f);
    segments = (seg *)malloc(4*nx*ny*sizeof(seg));
    if (!segments) return GEOMLIB_MESH_ERROR_MEMORY;

    for (y = 0; y + 1 < ny; ++y) {
        for (x = 0; x + 1 < nx; ++x) {
            size_t i00 = x + nx*y;
            size_t i10 = i00+1;
            size_t i01 = i00+nx;
            size_t i11 = i01+1;
            double v[4] = {values[i00], values[i10], values[i11], values[i01]};
            geomlib_vec2 p[4] = {
                {origin.x + (float)x*h, origin.y + (float)y*h},
                {origin.x + (float)(x+1)*h, origin.y + (float)y*h},
                {origin.x + (float)(x+1)*h, origin.y + (float)(y+1)*h},
                {origin.x + (float)x*h, origin.y + (float)(y+1)*h}
            };
            int inside[4], code = 0, e, n = 0;
            geomlib_vec2 list[4];
            for (e = 0; e < 4; ++e) {
                inside[e] = v[e] < 0.0;
                if (inside[e]) code |= (1 << e);
            }
            if (code == 0 || code == 15) continue;
            for (e = 0; e < 4; ++e) {
                int f = (e+1)%4;
                if (inside[e] == inside[f]) continue;
                list[n++] = lerp2(p[e], p[f], v[e], v[f]);
            }
            if (n == 2) {
                segments[nseg].a = list[0];
                segments[nseg].b = list[1];
                nseg++;
            } else if (n == 4) {
                double center = (v[0]+v[1]+v[2]+v[3])*0.25;
                int c = center < 0.0;
                if (c == inside[0]) {
                    segments[nseg].a=list[0]; segments[nseg++].b=list[1];
                    segments[nseg].a=list[2]; segments[nseg++].b=list[3];
                } else {
                    segments[nseg].a=list[1]; segments[nseg++].b=list[2];
                    segments[nseg].a=list[3]; segments[nseg++].b=list[0];
                }
            }
        }
    }
    if (!nseg) { free(segments); memset(out, 0, sizeof(*out)); return GEOMLIB_MESH_SUCCESS; }

    {
        int *used = (int *)calloc(nseg, sizeof(int));
        geomlib_vec2 *points = (geomlib_vec2 *)malloc(2*nseg*sizeof(geomlib_vec2));
        uint32_t *starts = (uint32_t *)malloc(nseg*sizeof(uint32_t));
        uint32_t *counts = (uint32_t *)calloc(nseg, sizeof(uint32_t));
        size_t np = 0, nlines = 0;
        if (!used || !points || !starts || !counts) {
            free(used); free(points); free(starts); free(counts); free(segments);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        for (i = 0; i < nseg; ++i) {
            geomlib_vec3 cur_a; /* unused */
            (void)cur_a;
            geomlib_vec2 cur, prev;
            size_t j;
            if (used[i]) continue;
            used[i] = 1;
            cur = segments[i].b;
            prev = segments[i].a;
            starts[nlines] = (uint32_t)np;
            points[np++] = prev;
            points[np++] = cur;
            counts[nlines] = 2;
            /* Extend at the tail. */
            for (;;) {
                int found = -1, at_b = 0;
                for (j = 0; j < nseg; ++j) {
                    if (used[j]) continue;
                    if (same(segments[j].a, cur, eps)) { found=(int)j; at_b=0; break; }
                    if (same(segments[j].b, cur, eps)) { found=(int)j; at_b=1; break; }
                }
                if (found < 0) break;
                used[found] = 1;
                cur = at_b ? segments[found].a : segments[found].b;
                if (!same(cur, points[np-1], eps)) {
                    points[np++] = cur;
                    counts[nlines]++;
                }
            }
            /* Prepend at the head. */
            for (;;) {
                int found = -1, at_b = 0;
                for (j = 0; j < nseg; ++j) {
                    if (used[j]) continue;
                    if (same(segments[j].b, prev, eps)) { found=(int)j; at_b=0; break; }
                    if (same(segments[j].a, prev, eps)) { found=(int)j; at_b=1; break; }
                }
                if (found < 0) break;
                used[found] = 1;
                prev = at_b ? segments[found].b : segments[found].a;
                {
                    size_t pos = starts[nlines];
                    size_t cnt = counts[nlines];
                    for (j = 0; j < cnt; ++j) points[pos+cnt-j] = points[pos+cnt-j-1];
                    points[pos] = prev;
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
        free(used);
    }
    free(segments);
    return GEOMLIB_MESH_SUCCESS;
}
