/*
 * geomlib - Ellipse queries (see geomlib_ellipse.h)
 *
 * Clean-room rotation into the ellipse frame. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_ellipse.h"

#include <math.h>

double geomlib_ellipse_implicit(geomlib_vec2 q, geomlib_vec2 center,
                                float a, float b, float angle)
{
    double dx = (double)q.x - center.x;
    double dy = (double)q.y - center.y;
    double c = cos(angle), s = sin(angle);
    double x = c * dx + s * dy;   /* rotate by -angle */
    double y = -s * dx + c * dy;
    double aa = fabs((double)a), bb = fabs((double)b);
    if (!(aa > 1e-12 && bb > 1e-12)) return INFINITY;
    return (x * x) / (aa * aa) + (y * y) / (bb * bb);
}

int geomlib_point_in_ellipse(geomlib_vec2 q, geomlib_vec2 center,
                             float a, float b, float angle)
{
    return geomlib_ellipse_implicit(q, center, a, b, angle) <= 1.0 + 1e-12;
}
