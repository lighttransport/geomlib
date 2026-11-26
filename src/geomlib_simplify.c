/*
 * geomlib - Polyline simplification (see geomlib_simplify.h)
 *
 * Clean-room Douglas-Peucker recursion. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_simplify.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static float point_segment_distance(geomlib_vec3 p, geomlib_vec3 a, geomlib_vec3 b)
{
    geomlib_vec3 ab = geomlib_vec3_sub(b,a);
    double len2 = geomlib_vec3_length_squared(ab);
    double t;
    geomlib_vec3 c;
    if (!(len2 > 1e-18)) return geomlib_vec3_length(geomlib_vec3_sub(p,a));
    t = geomlib_vec3_dot(geomlib_vec3_sub(p,a), ab)/len2;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    c = geomlib_vec3_add(a, geomlib_vec3_mul(ab, (float)t));
    return geomlib_vec3_length(geomlib_vec3_sub(p,c));
}

static void dp_recurse(const geomlib_vec3 *p, size_t first, size_t last,
                       float eps, uint8_t *keep)
{
    size_t i, best=0;
    float maxd = -1.0f;
    if (last <= first+1) return;
    for (i = first+1; i < last; ++i) {
        float d = point_segment_distance(p[i], p[first], p[last]);
        if (d > maxd) { maxd = d; best = i; }
    }
    if (maxd > eps) {
        keep[best] = 1;
        dp_recurse(p, first, best, eps, keep);
        dp_recurse(p, best, last, eps, keep);
    }
}

int geomlib_simplify_polyline(const geomlib_vec3 *points, size_t count,
                              float epsilon, uint32_t *indices,
                              size_t *out_count)
{
    uint8_t *keep;
    size_t i, n = 0;
    if (!points || !indices || !out_count || !count || !(epsilon >= 0.0f))
        return GEOMLIB_MESH_ERROR_INVALID;
    keep = (uint8_t *)calloc(count, sizeof(uint8_t));
    if (!keep) return GEOMLIB_MESH_ERROR_MEMORY;
    keep[0] = 1;
    if (count > 1) keep[count-1] = 1;
    dp_recurse(points, 0, count-1, epsilon, keep);
    for (i = 0; i < count; ++i)
        if (keep[i]) indices[n++] = (uint32_t)i;
    free(keep);
    *out_count = n;
    return GEOMLIB_MESH_SUCCESS;
}

float geomlib_simplify_max_error(const geomlib_vec3 *points, size_t count,
                                 const uint32_t *indices, size_t kept)
{
    size_t i, k;
    float mx = 0.0f;
    if (!points || !indices || !kept) return INFINITY;
    for (i = 0; i < count; ++i) {
        float best = INFINITY;
        for (k = 0; k+1 < kept; ++k) {
            float d = point_segment_distance(points[i],
                points[indices[k]], points[indices[k+1]]);
            if (d < best) best = d;
        }
        if (best > mx) mx = best;
    }
    return mx;
}
