/*
 * geomlib - Plane / sphere intersection (see geomlib_plane_sphere.h)
 *
 * Clean-room analytic circle computation. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_plane_sphere.h"

#include <float.h>
#include <math.h>
#include <string.h>

geomlib_plane_sphere_intersection geomlib_plane_sphere_intersect(
    geomlib_vec3 plane_point, geomlib_vec3 plane_normal,
    geomlib_vec3 sphere_center, float radius)
{
    geomlib_plane_sphere_intersection r;
    double n[3], len, d, h2;

    memset(&r, 0, sizeof(r));
    r.kind = GEOMLIB_PLANE_SPHERE_NONE;
    r.center = geomlib_vec3_create(0, 0, 0);
    r.normal = geomlib_vec3_create(0, 0, 0);
    r.radius = 0.0f;
    if (!(radius >= 0.0f)) return r;
    n[0] = (double)plane_normal.x;
    n[1] = (double)plane_normal.y;
    n[2] = (double)plane_normal.z;
    len = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (!(len > 1e-12)) return r;
    n[0] /= len; n[1] /= len; n[2] /= len;
    r.normal = geomlib_vec3_create((float)n[0], (float)n[1], (float)n[2]);

    d = (double)sphere_center.x * n[0] + (double)sphere_center.y * n[1]
      + (double)sphere_center.z * n[2]
      - ((double)plane_point.x * n[0] + (double)plane_point.y * n[1]
         + (double)plane_point.z * n[2]);
    if (fabs(d) > (double)radius + 1e-9) return r;
    r.center = geomlib_vec3_create(
        (float)((double)sphere_center.x - d * n[0]),
        (float)((double)sphere_center.y - d * n[1]),
        (float)((double)sphere_center.z - d * n[2]));
    h2 = (double)radius * radius - d * d;
    if (h2 < -1e-9) return r;
    r.radius = (float)sqrt(fmax(0.0, h2));
    r.kind = (r.radius < 1e-5f) ? GEOMLIB_PLANE_SPHERE_TANGENT
                                : GEOMLIB_PLANE_SPHERE_CIRCLE;
    return r;
}
