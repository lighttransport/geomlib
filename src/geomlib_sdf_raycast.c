/*
 * geomlib - SDF ray casting (see geomlib_sdf_raycast.h)
 *
 * Clean-room sphere tracing with sign-change bisection and central-difference
 * normals. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_sdf_raycast.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static float sdf_at(const geomlib_sdf *sdf, geomlib_vec3 p)
{
    return geomlib_sdf_query_distance(sdf, p);
}

int geomlib_sdf_raycast(const geomlib_sdf *sdf, geomlib_ray ray,
                        float max_distance, uint32_t max_steps,
                        geomlib_sdf_ray_hit *out)
{
    float t, eps, d;
    uint32_t step;
    geomlib_vec3 dir;

    if (!sdf || !out || !sdf->values) return GEOMLIB_MESH_ERROR_INVALID;
    memset(out, 0, sizeof(*out));
    if (!isfinite(ray.direction.x) || !isfinite(ray.direction.y) ||
        !isfinite(ray.direction.z)) return GEOMLIB_MESH_ERROR_INVALID;
    dir = geomlib_vec3_normalize(ray.direction);
    if (geomlib_vec3_length(dir) < 0.5f) return GEOMLIB_MESH_ERROR_INVALID;
    eps = sdf->voxel_size * 0.25f;
    if (!(eps > 1e-9f)) eps = 1e-6f;
    if (!(max_distance > 0.0f)) max_distance = sdf->extent.x +
        sdf->extent.y + sdf->extent.z;
    if (!max_steps) max_steps = 256;

    t = ray.t_min > 0.0f ? ray.t_min : 0.0f;
    for (step = 0; step < max_steps; ++step) {
        geomlib_vec3 p = geomlib_vec3_add(ray.origin, geomlib_vec3_mul(dir,t));
        d = sdf_at(sdf, p);
        out->steps = step+1;
        if (!isfinite(d)) break;
        if (fabs(d) <= eps) {
            out->hit = 1;
            break;
        }
        t += fmaxf(fabs(d), eps);
        if (t > max_distance) break;
    }
    if (!out->hit) return GEOMLIB_MESH_SUCCESS;

    out->t = t;
    out->position = geomlib_vec3_add(ray.origin, geomlib_vec3_mul(dir,t));
    {
        float h = sdf->voxel_size;
        float dx = sdf_at(sdf, geomlib_vec3_add(out->position,
                          geomlib_vec3_create(h,0,0)))
                 - sdf_at(sdf, geomlib_vec3_sub(out->position,
                          geomlib_vec3_create(h,0,0)));
        float dy = sdf_at(sdf, geomlib_vec3_add(out->position,
                          geomlib_vec3_create(0,h,0)))
                 - sdf_at(sdf, geomlib_vec3_sub(out->position,
                          geomlib_vec3_create(0,h,0)));
        float dz = sdf_at(sdf, geomlib_vec3_add(out->position,
                          geomlib_vec3_create(0,0,h)))
                 - sdf_at(sdf, geomlib_vec3_sub(out->position,
                          geomlib_vec3_create(0,0,h)));
        geomlib_vec3 n = geomlib_vec3_create(dx,dy,dz);
        float len = geomlib_vec3_length(n);
        out->normal = len > 1e-12f ? geomlib_vec3_mul(n,1.0f/len)
                                   : geomlib_vec3_create(0,0,1);
    }
    return GEOMLIB_MESH_SUCCESS;
}
