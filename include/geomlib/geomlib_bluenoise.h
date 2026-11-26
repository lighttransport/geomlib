/*
 * geomlib - Blue-noise surface sampling
 *
 * Bridson's Poisson-disk sampling generalized to triangle meshes: candidate
 * points are generated in triangles with area-weighted selection and rejected
 * if they are too close to an accepted sample (spatial hash acceleration).
 */

#ifndef GEOMLIB_BLUENOISE_H
#define GEOMLIB_BLUENOISE_H

#include "geomlib_mesh.h"
#include "geomlib_random.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 *points;
    geomlib_vec3 *normals;
    uint32_t *triangle_ids;
    size_t count;
} geomlib_sample_set;

geomlib_sample_set geomlib_sample_set_create(size_t capacity);
void geomlib_sample_set_free(geomlib_sample_set *set);

/*
 * Generate blue-noise samples on a triangle mesh. `min_distance` is the
 * minimum world-space separation; `rng` may be NULL (deterministic seed).
 * Returns the number of samples generated.
 */
size_t geomlib_blue_noise_surface(const geomlib_trimesh *mesh,
                                  float min_distance,
                                  geomlib_pcg32 *rng,
                                  geomlib_sample_set *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_BLUENOISE_H */
