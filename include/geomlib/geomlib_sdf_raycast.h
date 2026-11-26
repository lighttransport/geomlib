/*
 * geomlib - Ray casting against a signed distance field
 *
 * Sphere tracing with bisection refinement and gradient-based normals.
 */

#ifndef GEOMLIB_SDF_RAYCAST_H
#define GEOMLIB_SDF_RAYCAST_H

#include "geomlib_sdf.h"
#include "geomlib_bvh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int hit;                  /* 1 when a zero crossing was found */
    float t;                  /* Ray parameter of the surface point */
    geomlib_vec3 position;
    geomlib_vec3 normal;      /* SDF gradient (outward for negative-inside SDFs) */
    uint32_t steps;
} geomlib_sdf_ray_hit;

/* Cast `ray` through `sdf`. Returns GEOMLIB_MESH_SUCCESS when the query itself
 * succeeds; inspect `out.hit` for whether a surface was reached. */
int geomlib_sdf_raycast(const geomlib_sdf *sdf, geomlib_ray ray,
                        float max_distance, uint32_t max_steps,
                        geomlib_sdf_ray_hit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SDF_RAYCAST_H */
