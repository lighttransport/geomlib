/*
 * geomlib - Convex polygon diameter (see geomlib_polygon_diameter.h)
 *
 * Clean-room rotating calipers. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polygon_diameter.h"

#include <math.h>

static double cross2(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c)
{
    return ((double)b.x - a.x) * ((double)c.y - a.y)
         - ((double)b.y - a.y) * ((double)c.x - a.x);
}

static double dist2(geomlib_vec2 a, geomlib_vec2 b)
{
    double dx = (double)a.x - b.x, dy = (double)a.y - b.y;
    return dx * dx + dy * dy;
}

double geomlib_convex_polygon_diameter(const geomlib_vec2 *poly, size_t n)
{
    size_t i, j = 1;
    double best = 0.0;

    if (!poly || n < 2) return 0.0;
    if (n == 2) return sqrt(dist2(poly[0], poly[1]));

    for (i = 0; i < n; ++i) {
        size_t ni = (i + 1) % n;
        double cur = fabs(cross2(poly[i], poly[ni], poly[j]));
        for (;;) {
            size_t nj = (j + 1) % n;
            double nxt = fabs(cross2(poly[i], poly[ni], poly[nj]));
            if (nxt <= cur + 1e-12) break;
            j = nj;
            cur = nxt;
            if (j == i) break;
        }
        {
            double d = dist2(poly[i], poly[j]);
            if (d > best) best = d;
            d = dist2(poly[ni], poly[j]);
            if (d > best) best = d;
        }
    }
    return sqrt(best);
}
