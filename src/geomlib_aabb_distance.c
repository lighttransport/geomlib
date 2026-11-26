/*
 * geomlib - AABB / AABB distance (see geomlib_aabb_distance.h)
 *
 * Clean-room per-axis gap decomposition. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_aabb_distance.h"

#include <math.h>

int geomlib_aabb_distance(geomlib_vec3 amin, geomlib_vec3 amax,
                          geomlib_vec3 bmin, geomlib_vec3 bmax,
                          geomlib_aabb_distance_result *out)
{
    geomlib_aabb_distance_result r;
    double pa[3], pb[3], sum = 0.0;
    double a0[3] = {(double)amin.x, (double)amin.y, (double)amin.z};
    double a1[3] = {(double)amax.x, (double)amax.y, (double)amax.z};
    double b0[3] = {(double)bmin.x, (double)bmin.y, (double)bmin.z};
    double b1[3] = {(double)bmax.x, (double)bmax.y, (double)bmax.z};
    int i;

    r.distance = 0.0f;
    r.point_a = geomlib_vec3_create(0, 0, 0);
    r.point_b = geomlib_vec3_create(0, 0, 0);
    if (!out) return 0;
    for (i = 0; i < 3; ++i) {
        double lo = fmax(a0[i], b0[i]);
        double hi = fmin(a1[i], b1[i]);
        double gap;
        if (a1[i] < b0[i]) {
            gap = b0[i] - a1[i];
            pa[i] = a1[i];
            pb[i] = b0[i];
        } else if (b1[i] < a0[i]) {
            gap = a0[i] - b1[i];
            pa[i] = a0[i];
            pb[i] = b1[i];
        } else {
            gap = 0.0;
            pa[i] = pb[i] = 0.5 * (lo + hi);
        }
        sum += gap * gap;
    }
    r.point_a = geomlib_vec3_create((float)pa[0], (float)pa[1], (float)pa[2]);
    r.point_b = geomlib_vec3_create((float)pb[0], (float)pb[1], (float)pb[2]);
    r.distance = (float)sqrt(sum);
    *out = r;
    return 1;
}
