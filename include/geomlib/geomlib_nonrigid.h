/*
 * geomlib - Non-rigid ICP registration
 *
 * Iteratively matches source vertices to closest points on a target surface
 * and fits a smoothed per-vertex displacement field (embedded deformation).
 */

#ifndef GEOMLIB_NONRIGID_H
#define GEOMLIB_NONRIGID_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 *displacement;   /* one offset per source vertex */
    size_t vertex_count;
    size_t iterations;
    float initial_rms;
    float final_rms;
} geomlib_nonrigid_result;

geomlib_nonrigid_result geomlib_nonrigid_result_create(void);
void geomlib_nonrigid_result_free(geomlib_nonrigid_result *r);

/* Fit a smooth deformation of `source` onto `target`.
 * `smoothness` in [0,1] blends per-vertex matches with neighbour averages
 * (0 = fully smooth, 1 = no smoothing). */
int geomlib_nonrigid_icp(const geomlib_trimesh *source,
                         const geomlib_trimesh *target,
                         uint32_t max_iterations, float smoothness,
                         geomlib_nonrigid_result *out);

/* Apply a fitted deformation, producing a new owned mesh. */
int geomlib_nonrigid_apply(const geomlib_trimesh *source,
                           const geomlib_nonrigid_result *fit,
                           geomlib_trimesh *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_NONRIGID_H */
