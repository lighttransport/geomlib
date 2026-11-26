/*
 * geomlib - Circle / line intersection (see geomlib_circle_line.h)
 *
 * Clean-room analytic quadratic. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_circle_line.h"

#include <math.h>
#include <string.h>

static geomlib_circle_line_intersection empty(void)
{
    geomlib_circle_line_intersection r;
    r.kind = GEOMLIB_CIRCLE_LINE_NONE;
    r.p0 = r.p1 = geomlib_vec2_create(0, 0);
    return r;
}

static geomlib_circle_line_intersection solve(geomlib_vec2 a, geomlib_vec2 b,
                                              geomlib_vec2 center, float radius,
                                              int finite)
{
    geomlib_circle_line_intersection r = empty();
    double ux = (double)b.x - a.x, uy = (double)b.y - a.y;
    double vx = (double)a.x - center.x, vy = (double)a.y - center.y;
    double aq = ux * ux + uy * uy, bq, cq, disc, sq, t0, t1;
    geomlib_vec2 pts[2];
    int np = 0;

    if (!(radius >= 0.0f) || !(aq > 1e-18)) return r;
    bq = vx * ux + vy * uy;
    cq = vx * vx + vy * vy - (double)radius * radius;
    disc = bq * bq - aq * cq;
    if (disc < 0.0) return r;
    sq = sqrt(disc);
    t0 = (-bq - sq) / aq;
    t1 = (-bq + sq) / aq;
    {
        double ts[2] = {t0, t1};
        int i;
        for (i = 0; i < 2; ++i) {
            if (finite && (ts[i] < -1e-12 || ts[i] > 1.0 + 1e-12)) continue;
            pts[np++] = geomlib_vec2_create(
                (float)((double)a.x + ts[i] * ux),
                (float)((double)a.y + ts[i] * uy));
        }
    }
    if (!np) return r;
    r.p0 = pts[0];
    r.p1 = pts[np - 1];
    r.kind = (np == 1) ? GEOMLIB_CIRCLE_LINE_TANGENT
                       : GEOMLIB_CIRCLE_LINE_TWO_POINTS;
    if (np == 2 && fabs(r.p0.x - r.p1.x) < 1e-5 &&
        fabs(r.p0.y - r.p1.y) < 1e-5) {
        r.kind = GEOMLIB_CIRCLE_LINE_TANGENT;
        r.p1 = r.p0;
    }
    return r;
}

geomlib_circle_line_intersection geomlib_circle_line_intersect(
    geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 center, float radius)
{
    return solve(a, b, center, radius, 0);
}

geomlib_circle_line_intersection geomlib_circle_segment_intersect(
    geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 center, float radius)
{
    return solve(a, b, center, radius, 1);
}
