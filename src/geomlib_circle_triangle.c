/*
 * geomlib - Circle / triangle intersection (see geomlib_circle_triangle.h)
 *
 * Clean-room containment and edge-distance test. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_circle_triangle.h"
#include "geomlib/geomlib_point_in_polygon.h"

#include <float.h>
#include <math.h>

static void closest_on_segment(geomlib_vec2 q, geomlib_vec2 a, geomlib_vec2 b,
                               geomlib_vec2 *out)
{
    double ux = (double)b.x - a.x, uy = (double)b.y - a.y;
    double u2 = ux * ux + uy * uy;
    double t = 0.0;
    if (u2 > 1e-18) {
        t = (((double)q.x - a.x) * ux + ((double)q.y - a.y) * uy) / u2;
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
    }
    out->x = (float)((double)a.x + t * ux);
    out->y = (float)((double)a.y + t * uy);
}

int geomlib_circle_triangle_intersect(geomlib_vec2 center, float radius,
                                      geomlib_vec2 a, geomlib_vec2 b,
                                      geomlib_vec2 c,
                                      geomlib_circle_triangle_result *out)
{
    geomlib_circle_triangle_result r;
    geomlib_vec2 tri[3], q;
    double best = DBL_MAX, dx, dy, d2;
    int i, all_inside = 1, any_vertex = 0;

    r.intersects = 0;
    r.center_inside = 0;
    r.contains_triangle = 0;
    r.signed_distance = 0.0f;
    r.closest = geomlib_vec2_create(0, 0);
    if (!out || !(radius >= 0.0f)) return 0;
    tri[0] = a; tri[1] = b; tri[2] = c;
    r.center_inside = geomlib_point_in_polygon(tri, 3, center) ? 1 : 0;
    for (i = 0; i < 3; ++i) {
        dx = (double)tri[i].x - center.x;
        dy = (double)tri[i].y - center.y;
        d2 = dx * dx + dy * dy;
        if (d2 <= (double)radius * radius + 1e-12) any_vertex = 1;
        else all_inside = 0;
    }
    for (i = 0; i < 3; ++i) {
        closest_on_segment(center, tri[i], tri[(i + 1) % 3], &q);
        dx = (double)center.x - q.x;
        dy = (double)center.y - q.y;
        d2 = dx * dx + dy * dy;
        if (d2 < best) { best = d2; r.closest = q; }
    }
    r.signed_distance = (float)(sqrt(best) - radius);
    r.contains_triangle = all_inside;
    r.intersects = (r.center_inside || any_vertex || r.signed_distance <= 0.0f);
    *out = r;
    return r.intersects;
}
