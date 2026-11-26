/*
 * geomlib - Variational Shape Approximation (VSA)
 *
 * Clusters triangles into planar proxy regions by minimizing a sum of
 * squared plane distances, then emits one fan patch per region boundary loop.
 */

#ifndef GEOMLIB_VSA_H
#define GEOMLIB_VSA_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t *region_of_triangle;  /* Cluster id per input triangle */
    size_t triangle_count;
    geomlib_vec3 *proxy_points;    /* Representative point per cluster */
    geomlib_vec3 *proxy_normals;   /* Fitted plane normal per cluster */
    size_t region_count;
} geomlib_vsa_result;

geomlib_vsa_result geomlib_vsa_result_create(void);
void geomlib_vsa_result_free(geomlib_vsa_result *r);

/*
 * Approximate `input` with `region_count` planar proxies and return both the
 * clustering and a simplified triangle mesh built from the region fans.
 */
int geomlib_vsa_simplify(const geomlib_trimesh *input,
                         size_t region_count,
                         int iterations,
                         geomlib_vsa_result *result,
                         geomlib_trimesh *output);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_VSA_H */
