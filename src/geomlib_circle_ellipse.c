/*
 * geomlib - Circle / ellipse overlap (see geomlib_circle_ellipse.h)
 *
 * Clean-room closest-point test. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_circle_ellipse.h"
#include "geomlib/geomlib_ellipse.h"
#include "geomlib/geomlib_ellipse_closest.h"

int geomlib_circle_ellipse_overlap(geomlib_vec2 circle_center, float radius,
                                   geomlib_vec2 ellipse_center,
                                   float a, float b, float angle,
                                   geomlib_circle_ellipse_result *out)
{
    geomlib_circle_ellipse_result r;
    geomlib_ellipse_closest cp;

    r.intersects = 0;
    r.center_inside = 0;
    r.distance = 0.0f;
    r.closest = geomlib_vec2_create(0, 0);
    if (!out || !(radius >= 0.0f)) return 0;
    if (!geomlib_ellipse_closest_point(circle_center, ellipse_center,
                                       a, b, angle, &cp)) return 0;
    r.center_inside = geomlib_point_in_ellipse(
        circle_center, ellipse_center, a, b, angle) ? 1 : 0;
    r.closest = cp.point;
    r.distance = cp.distance - radius;
    r.intersects = (r.center_inside || r.distance <= 0.0f);
    *out = r;
    return r.intersects;
}
