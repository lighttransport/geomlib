/*
 * geomlib - Ray / OBB intersection (see geomlib_ray_obb.h)
 *
 * Clean-room slab test in local coordinates. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_ray_obb.h"

#include <float.h>
#include <math.h>

int geomlib_ray_obb_intersect(geomlib_vec3 origin, geomlib_vec3 direction,
                              const geomlib_obb *obb,
                              geomlib_ray_obb_hit *out)
{
    double o[3], d[3], lo[3], ld[3], h[3];
    double tmin = -DBL_MAX, tmax = DBL_MAX;
    int near_axis = 0, near_sign = 1, far_axis = 0, far_sign = 1, i;
    geomlib_ray_obb_hit r;

    r.t = 0.0;
    r.point = geomlib_vec3_create(0, 0, 0);
    r.normal = geomlib_vec3_create(0, 0, 0);
    if (!out || !obb) return 0;
    o[0] = (double)origin.x; o[1] = (double)origin.y; o[2] = (double)origin.z;
    d[0] = (double)direction.x; d[1] = (double)direction.y;
    d[2] = (double)direction.z;
    for (i = 0; i < 3; ++i) {
        geomlib_vec3 ax = obb->axes[i];
        double axn = sqrt((double)ax.x * ax.x + (double)ax.y * ax.y +
                          (double)ax.z * ax.z);
        if (!(axn > 1e-12)) return 0;
        lo[i] = ((double)origin.x - obb->center.x) * ax.x +
                ((double)origin.y - obb->center.y) * ax.y +
                ((double)origin.z - obb->center.z) * ax.z;
        ld[i] = d[0] * ax.x + d[1] * ax.y + d[2] * ax.z;
        h[i] = fabs(i == 0 ? (double)obb->half_extents.x
                    : (i == 1 ? (double)obb->half_extents.y
                              : (double)obb->half_extents.z));
    }
    for (i = 0; i < 3; ++i) {
        double t0, t1, sign;
        if (fabs(ld[i]) < 1e-15) {
            if (lo[i] < -h[i] || lo[i] > h[i]) return 0;
            continue;
        }
        t0 = (-h[i] - lo[i]) / ld[i];
        t1 = (h[i] - lo[i]) / ld[i];
        sign = 1;
        if (t0 > t1) { double t = t0; t0 = t1; t1 = t; sign = -1; }
        if (t0 > tmin) { tmin = t0; near_axis = i; near_sign = sign; }
        if (t1 < tmax) { tmax = t1; far_axis = i; far_sign = sign; }
        if (tmin > tmax) return 0;
    }
    if (tmax < 0.0) return 0;
    {
        double t = tmin >= 0.0 ? tmin : tmax;
        int axis = tmin >= 0.0 ? near_axis : far_axis;
        int sign = tmin >= 0.0 ? -near_sign : far_sign;
        geomlib_vec3 ax = obb->axes[axis];
        r.t = t;
        r.point = geomlib_vec3_create(
            (float)(o[0] + t * d[0]),
            (float)(o[1] + t * d[1]),
            (float)(o[2] + t * d[2]));
        r.normal = geomlib_vec3_mul(ax, (float)sign);
    }
    *out = r;
    return 1;
}
