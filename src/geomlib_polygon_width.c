/*
 * geomlib - Convex polygon minimum width (see geomlib_polygon_width.h)
 *
 * Clean-room rotating calipers on antipodal vertex-edge pairs. C11, no
 * external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polygon_width.h"

#include <float.h>
#include <math.h>

static double cross2(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 c)
{
    return ((double)b.x - a.x) * ((double)c.y - a.y)
         - ((double)b.y - a.y) * ((double)c.x - a.x);
}

double geomlib_convex_polygon_min_width(const geomlib_vec2 *poly, size_t n)
{
    size_t i, j = 1;
    double best = DBL_MAX;

    if (!poly || n < 3) return 0.0;
    for (i = 0; i < n; ++i) {
        size_t ni = (i + 1) % n;
        double ex = (double)poly[ni].x - poly[i].x;
        double ey = (double)poly[ni].y - poly[i].y;
        double elen = sqrt(ex * ex + ey * ey);
        double cur;
        if (!(elen > 1e-12)) continue;
        cur = fabs(cross2(poly[i], poly[ni], poly[j]));
        for (;;) {
            size_t nj = (j + 1) % n;
            double nxt = fabs(cross2(poly[i], poly[ni], poly[nj]));
            if (nxt <= cur + 1e-12) break;
            j = nj;
            cur = nxt;
            if (j == i) break;
        }
        {
            double w = cur / elen;
            if (w < best) best = w;
        }
    }
    return best == DBL_MAX ? 0.0 : best;
}
