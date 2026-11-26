/*
 * geomlib - Line / line intersection (see geomlib_line2d.h)
 *
 * Clean-room cross-product solution. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_line2d.h"

#include <math.h>
#include <string.h>

static double cross(geomlib_vec2 u, geomlib_vec2 v)
{
    return (double)u.x * v.y - (double)u.y * v.x;
}

geomlib_line2d_intersection geomlib_line2d_intersect(
    geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c, geomlib_vec2 d)
{
    geomlib_line2d_intersection r;
    geomlib_vec2 u = geomlib_vec2_sub(b, a), v = geomlib_vec2_sub(d, c);
    geomlib_vec2 w = geomlib_vec2_sub(c, a);
    double den = cross(u, v), t, q;

    memset(&r, 0, sizeof(r));
    r.kind = GEOMLIB_LINE2D_NONE;
    r.point = geomlib_vec2_create(0, 0);
    if (fabs(den) < 1e-14) {
        if (fabs(cross(w, u)) < 1e-12) r.kind = GEOMLIB_LINE2D_COLLINEAR;
        return r;
    }
    t = cross(w, v) / den;
    q = cross(w, u) / den;
    r.kind = GEOMLIB_LINE2D_POINT;
    r.t = t;
    r.u = q;
    r.point = geomlib_vec2_create(
        (float)((double)a.x + t * u.x), (float)((double)a.y + t * u.y));
    return r;
}
