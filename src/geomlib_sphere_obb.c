/*
 * geomlib - Sphere / OBB overlap (see geomlib_sphere_obb.h)
 *
 * Clean-room local-frame clamp. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_sphere_obb.h"

#include <math.h>

static double comp(geomlib_vec3 v, int i)
{
    return i == 0 ? (double)v.x : (i == 1 ? (double)v.y : (double)v.z);
}

int geomlib_sphere_obb_overlap(geomlib_vec3 center, float radius,
                               const geomlib_obb *obb,
                               geomlib_sphere_obb_result *out)
{
    geomlib_sphere_obb_result r;
    double local[3], h[3], clamped[3];
    geomlib_vec3 closest;
    double dx, dy, dz, d;
    int i;

    r.intersects = 0;
    r.distance = 0.0f;
    r.closest = geomlib_vec3_create(0, 0, 0);
    if (!out || !obb || !(radius >= 0.0f)) return 0;
    for (i = 0; i < 3; ++i) {
        geomlib_vec3 ax = obb->axes[i];
        double axn = sqrt((double)ax.x * ax.x + (double)ax.y * ax.y +
                          (double)ax.z * ax.z);
        if (!(axn > 1e-12)) return 0;
        local[i] = ((double)center.x - obb->center.x) * ax.x +
                   ((double)center.y - obb->center.y) * ax.y +
                   ((double)center.z - obb->center.z) * ax.z;
        h[i] = fabs(comp(obb->half_extents, i));
        clamped[i] = local[i] < -h[i] ? -h[i] : (local[i] > h[i] ? h[i] : local[i]);
    }
    closest = obb->center;
    for (i = 0; i < 3; ++i)
        closest = geomlib_vec3_add(closest, geomlib_vec3_mul(obb->axes[i],
                                   (float)clamped[i]));
    dx = (double)center.x - closest.x;
    dy = (double)center.y - closest.y;
    dz = (double)center.z - closest.z;
    d = sqrt(dx * dx + dy * dy + dz * dz);
    r.closest = closest;
    r.distance = (float)(d - radius);
    r.intersects = (r.distance <= 0.0f);
    *out = r;
    return r.intersects;
}
