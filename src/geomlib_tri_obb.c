/*
 * geomlib - Triangle / OBB overlap (see geomlib_tri_obb.h)
 *
 * Clean-room local-frame transform plus the shared triangle-AABB SAT. C11,
 * no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_tri_obb.h"
#include "geomlib/geomlib_tri_aabb.h"

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

int geomlib_triangle_obb_overlap(geomlib_vec3 a0, geomlib_vec3 a1,
                                 geomlib_vec3 a2, const geomlib_obb *obb)
{
    geomlib_vec3 lo, hi;
    if (!obb) return 0;
    lo = geomlib_vec3_create(-fabsf(obb->half_extents.x),
                             -fabsf(obb->half_extents.y),
                             -fabsf(obb->half_extents.z));
    hi = geomlib_vec3_create(fabsf(obb->half_extents.x),
                             fabsf(obb->half_extents.y),
                             fabsf(obb->half_extents.z));
    return geomlib_triangle_aabb_overlap(
        to_local(a0, obb), to_local(a1, obb), to_local(a2, obb), lo, hi);
}
