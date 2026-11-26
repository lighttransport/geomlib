/*
 * geomlib - Circle / circle intersection (see geomlib_circle_circle.h)
 *
 * Clean-room analytic radical-line solution. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_circle_circle.h"

#include <float.h>
#include <math.h>
#include <string.h>

geomlib_circle2d_intersection geomlib_circle2d_intersect(
    geomlib_vec2 c0, float r0, geomlib_vec2 c1, float r1)
{
    geomlib_circle2d_intersection r;
    double dx, dy, d, a, h2, h, ux, uy, bx, by;
    const double eps = 1e-9;

    memset(&r, 0, sizeof(r));
    r.kind = GEOMLIB_CIRCLE2D_NONE;
    r.p0 = r.p1 = geomlib_vec2_create(0, 0);
    if (!(r0 >= 0.0f && r1 >= 0.0f)) return r;

    dx = (double)c1.x - c0.x;
    dy = (double)c1.y - c0.y;
    d = sqrt(dx * dx + dy * dy);
    if (d < eps) {
        if (fabs((double)r0 - r1) < eps) {
            r.kind = GEOMLIB_CIRCLE2D_COINCIDENT;
            r.p0 = r.p1 = c0;
        } else {
            r.kind = GEOMLIB_CIRCLE2D_CONTAINED;
        }
        return r;
    }
    ux = dx / d;
    uy = dy / d;
    if (d > (double)r0 + r1 + eps) return r;
    if (d < fabs((double)r0 - r1) - eps) {
        r.kind = GEOMLIB_CIRCLE2D_CONTAINED;
        return r;
    }

    a = ((double)r0 * r0 - (double)r1 * r1 + d * d) / (2.0 * d);
    h2 = (double)r0 * r0 - a * a;
    if (h2 < -eps) return r;
    h = sqrt(fmax(0.0, h2));
    bx = (double)c0.x + a * ux;
    by = (double)c0.y + a * uy;
    r.p0 = geomlib_vec2_create((float)(bx - h * uy), (float)(by + h * ux));
    r.p1 = geomlib_vec2_create((float)(bx + h * uy), (float)(by - h * ux));
    r.kind = (h < 1e-6) ? GEOMLIB_CIRCLE2D_TANGENT
                        : GEOMLIB_CIRCLE2D_TWO_POINTS;
    return r;
}
