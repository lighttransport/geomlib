/*
 * geomlib - Capsule / OBB overlap (see geomlib_capsule_obb.h)
 *
 * Clean-room local-frame transform plus shared segment-box distance. C11, no
 * external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_capsule_obb.h"
#include "geomlib/geomlib_capsule_aabb.h"

#include <math.h>

static geomlib_vec3 to_local(geomlib_vec3 p, const geomlib_obb *obb)
{
    double v[3] = {0.0, 0.0, 0.0};
    int i;
    for (i = 0; i < 3; ++i) {
        geomlib_vec3 ax = obb->axes[i];
        double axn = sqrt((double)ax.x * ax.x + (double)ax.y * ax.y +
                          (double)ax.z * ax.z);
        if (!(axn > 1e-12)) return geomlib_vec3_create(0, 0, 0);
        v[i] = ((double)p.x - obb->center.x) * ax.x +
               ((double)p.y - obb->center.y) * ax.y +
               ((double)p.z - obb->center.z) * ax.z;
    }
    return geomlib_vec3_create((float)v[0], (float)v[1], (float)v[2]);
}

int geomlib_capsule_obb_overlap(geomlib_vec3 a, geomlib_vec3 b, float r,
                                const geomlib_obb *obb,
                                geomlib_capsule_obb_result *out)
{
    geomlib_capsule_obb_result res;
    geomlib_vec3 la, lb, lo, hi;
    float d;

    res.intersects = 0;
    res.distance = 0.0f;
    if (!out || !obb || !(r >= 0.0f)) return 0;
    la = to_local(a, obb);
    lb = to_local(b, obb);
    lo = geomlib_vec3_create(-fabsf(obb->half_extents.x),
                             -fabsf(obb->half_extents.y),
                             -fabsf(obb->half_extents.z));
    hi = geomlib_vec3_create(fabsf(obb->half_extents.x),
                             fabsf(obb->half_extents.y),
                             fabsf(obb->half_extents.z));
    d = geomlib_segment_aabb_distance(la, lb, lo, hi);
    res.distance = d - r;
    res.intersects = (res.distance <= 0.0f);
    *out = res;
    return res.intersects;
}
