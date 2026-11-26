/*
 * geomlib - Ray / tetrahedron intersection (see geomlib_ray_tetra.h)
 *
 * Clean-room enumeration of the four face intersections using the shared
 * ray/triangle test. C11, wasm32 safe.
 */

#include "geomlib/geomlib_ray_tetra.h"
#include "geomlib/geomlib_bvh.h"

#include <float.h>

static void face_at(const geomlib_vec3 v[4], int f,
                    geomlib_vec3 *a, geomlib_vec3 *b, geomlib_vec3 *c)
{
    static const int idx[4][3] = {{0, 1, 2}, {0, 3, 1}, {0, 2, 3}, {1, 3, 2}};
    *a = v[idx[f][0]];
    *b = v[idx[f][1]];
    *c = v[idx[f][2]];
}

int geomlib_ray_tetrahedron_intersect(
    geomlib_vec3 origin, geomlib_vec3 direction,
    geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2, geomlib_vec3 v3,
    geomlib_ray_tetra_hit *out)
{
    geomlib_vec3 v[4], a, b, c;
    geomlib_ray ray;
    geomlib_ray_tetra_hit r;
    float best = FLT_MAX, dlen;
    int f, found = -1;

    r.intersects = 0;
    r.t = 0.0f;
    r.point = geomlib_vec3_create(0, 0, 0);
    r.normal = geomlib_vec3_create(0, 0, 0);
    if (!out) return 0;
    v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
    dlen = geomlib_vec3_length(direction);
    if (!(dlen > 1e-12f)) return 0;
    ray = geomlib_ray_create(origin, direction);
    for (f = 0; f < 4; ++f) {
        float t = 0.0f, u = 0.0f, w = 0.0f;
        face_at(v, f, &a, &b, &c);
        if (geomlib_ray_triangle_intersect(&ray, a, b, c, &t, &u, &w) &&
            t >= 0.0f && t / dlen < best) {
            best = t / dlen;
            found = f;
        }
    }
    if (found < 0) return 0;
    face_at(v, found, &a, &b, &c);
    r.intersects = 1;
    r.t = best;
    r.point = geomlib_vec3_add(origin, geomlib_vec3_mul(direction, best));
    {
        geomlib_vec3 n = geomlib_vec3_cross(
            geomlib_vec3_sub(b, a), geomlib_vec3_sub(c, a));
        float len = geomlib_vec3_length(n);
        r.normal = len > 1e-12f ? geomlib_vec3_mul(n, 1.0f / len)
                                : geomlib_vec3_create(0, 0, 0);
    }
    *out = r;
    return 1;
}
