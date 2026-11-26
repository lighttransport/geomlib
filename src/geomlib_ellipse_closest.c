/*
 * geomlib - Ellipse closest point (see geomlib_ellipse_closest.h)
 *
 * Clean-room multi-start Newton iteration. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_ellipse_closest.h"

#include <float.h>
#include <math.h>

static const double PI2 = 6.28318530717958647692;

static void eval(double t, double c, double s, double aa, double bb,
                 geomlib_vec2 center, double e[2], double ep[2])
{
    double u = aa * cos(t), v = bb * sin(t);
    double up = -aa * sin(t), vp = bb * cos(t);
    e[0] = center.x + c * u - s * v;
    e[1] = center.y + s * u + c * v;
    ep[0] = c * up - s * vp;
    ep[1] = s * up + c * vp;
}

int geomlib_ellipse_closest_point(geomlib_vec2 q, geomlib_vec2 center,
                                  float a, float b, float angle,
                                  geomlib_ellipse_closest *out)
{
    double aa = fabs((double)a), bb = fabs((double)b);
    double c = cos(angle), s = sin(angle);
    double best = DBL_MAX, best_t = 0.0, best_p[2] = {0, 0};
    int start, it;

    if (!out || !(aa > 1e-12 && bb > 1e-12)) return 0;
    for (start = 0; start < 8; ++start) {
        double t = PI2 * start / 8.0;
        for (it = 0; it < 40; ++it) {
            double e[2], ep[2], dx, dy, f, fp, step;
            eval(t, c, s, aa, bb, center, e, ep);
            dx = e[0] - q.x;
            dy = e[1] - q.y;
            f = dx * ep[0] + dy * ep[1];
            /* E'' = -(E - centre), so f' = |E'|^2 + (E-q).E'' */
            fp = ep[0] * ep[0] + ep[1] * ep[1]
               - (dx * (e[0] - center.x) + dy * (e[1] - center.y));
            if (fabs(fp) < 1e-18) break;
            step = f / fp;
            if (step > 0.3) step = 0.3;
            if (step < -0.3) step = -0.3;
            t -= step;
            if (fabs(step) < 1e-13) break;
        }
        {
            double e[2], ep[2], dx, dy, d2;
            eval(t, c, s, aa, bb, center, e, ep);
            dx = e[0] - q.x; dy = e[1] - q.y;
            d2 = dx * dx + dy * dy;
            if (d2 < best) {
                best = d2;
                best_t = t;
                best_p[0] = e[0];
                best_p[1] = e[1];
            }
        }
    }
    out->point = geomlib_vec2_create((float)best_p[0], (float)best_p[1]);
    out->distance = (float)sqrt(best);
    out->t = (float)fmod(fmod(best_t, PI2) + PI2, PI2);
    return 1;
}
