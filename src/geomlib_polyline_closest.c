/*
 * geomlib - Polyline closest point (see geomlib_polyline_closest.h)
 *
 * Clean-room per-segment projection. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polyline_closest.h"

#include <float.h>

static float clamp01(float v)
{
    if (v < 0.0f) return 0.0f;
    if (v > 1.0f) return 1.0f;
    return v;
}

int geomlib_polyline_closest_point(const geomlib_vec3 *points, size_t n,
                                   geomlib_vec3 q,
                                   geomlib_polyline_closest *out)
{
    geomlib_polyline_closest r;
    size_t i;

    r.point = geomlib_vec3_create(0, 0, 0);
    r.distance = FLT_MAX;
    r.segment = 0;
    r.t = 0.0f;
    if (!out || !points || !n) return 0;
    if (n == 1) {
        r.point = points[0];
        r.distance = geomlib_vec3_length(geomlib_vec3_sub(q, points[0]));
        *out = r;
        return 1;
    }
    for (i = 0; i + 1 < n; ++i) {
        geomlib_vec3 u = geomlib_vec3_sub(points[i + 1], points[i]);
        float u2 = geomlib_vec3_dot(u, u);
        float t = 0.0f;
        geomlib_vec3 p;
        float d;
        if (u2 > 1e-12f)
            t = clamp01(geomlib_vec3_dot(geomlib_vec3_sub(q, points[i]), u) / u2);
        p = geomlib_vec3_add(points[i], geomlib_vec3_mul(u, t));
        d = geomlib_vec3_length(geomlib_vec3_sub(q, p));
        if (d < r.distance) {
            r.distance = d;
            r.point = p;
            r.segment = i;
            r.t = t;
        }
    }
    *out = r;
    return 1;
}
