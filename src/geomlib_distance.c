/*
 * geomlib - Bounded-error Hausdorff distance (see geomlib_distance.h)
 *
 * Clean-room implementation from the 1-Lipschitz property of the unsigned
 * distance function to a closed set. C11, wasm32 compatible.
 */

#include "geomlib/geomlib_distance.h"
#include "geomlib/geomlib_bvh.h"

#include <math.h>

static float vec_dist(geomlib_vec3 a, geomlib_vec3 b)
{
    geomlib_vec3 d = geomlib_vec3_sub(a, b);
    return geomlib_vec3_length(d);
}

static void direction_bounds(const geomlib_trimesh *src,
                             const geomlib_bvh *target_bvh,
                             float *lower, float *upper, size_t *samples)
{
    size_t t;
    *lower = 0.0f;
    *upper = 0.0f;
    *samples = 0;

    for (t = 0; t < src->triangle_count; ++t) {
        geomlib_vec3 p[3], mid[3], centroid;
        float dv[3], dm[3], dc;
        float longest;
        float min_vertex, tri_upper;
        size_t k;

        for (k = 0; k < 3; ++k) {
            p[k] = src->vertices[src->indices[t * 3 + k]];
        }
        longest = vec_dist(p[0], p[1]);
        {
            float e = vec_dist(p[1], p[2]);
            if (e > longest) longest = e;
            e = vec_dist(p[2], p[0]);
            if (e > longest) longest = e;
        }
        for (k = 0; k < 3; ++k) {
            geomlib_vec3 closest;
            float distance;
            if (geomlib_bvh_closest_point(target_bvh, p[k], &closest, &distance, NULL)) {
                dv[k] = distance;
            } else {
                dv[k] = INFINITY;
            }
            mid[k] = geomlib_vec3_mul(geomlib_vec3_add(p[k], p[(k + 1) % 3]), 0.5f);
            if (geomlib_bvh_closest_point(target_bvh, mid[k], &closest, &distance, NULL)) {
                dm[k] = distance;
            } else {
                dm[k] = INFINITY;
            }
            *samples += 2;
        }
        centroid = geomlib_vec3_mul(
            geomlib_vec3_add(geomlib_vec3_add(p[0], p[1]), p[2]), 1.0f / 3.0f);
        {
            geomlib_vec3 closest;
            float distance;
            if (geomlib_bvh_closest_point(target_bvh, centroid, &closest, &distance, NULL)) {
                dc = distance;
            } else {
                dc = INFINITY;
            }
        }
        (*samples)++;

        min_vertex = dv[0] < dv[1] ? dv[0] : dv[1];
        if (dv[2] < min_vertex) min_vertex = dv[2];
        tri_upper = min_vertex + longest;

        for (k = 0; k < 3; ++k) {
            if (dv[k] > *lower) *lower = dv[k];
            if (dm[k] > *lower) *lower = dm[k];
        }
        if (dc > *lower) *lower = dc;
        if (tri_upper > *upper) *upper = tri_upper;
    }
}

int geomlib_hausdorff_bounded(const geomlib_trimesh *a,
                              const geomlib_trimesh *b,
                              geomlib_hausdorff_result *out)
{
    geomlib_bvh bvh_a, bvh_b;
    float lower_a, upper_a, lower_b, upper_b;
    size_t samples_a, samples_b;

    if (!a || !b || !out) return GEOMLIB_MESH_ERROR_INVALID;
    out->lower_bound = 0.0f;
    out->upper_bound = 0.0f;
    out->samples = 0;
    if (a->triangle_count == 0 && b->triangle_count == 0) return GEOMLIB_MESH_SUCCESS;
    if (a->triangle_count == 0 || b->triangle_count == 0) {
        out->lower_bound = INFINITY;
        out->upper_bound = INFINITY;
        return GEOMLIB_MESH_SUCCESS;
    }

    geomlib_bvh_build(&bvh_a, a, GEOMLIB_BVH_QUALITY_MEDIUM);
    geomlib_bvh_build(&bvh_b, b, GEOMLIB_BVH_QUALITY_MEDIUM);
    direction_bounds(a, &bvh_b, &lower_a, &upper_a, &samples_a);
    direction_bounds(b, &bvh_a, &lower_b, &upper_b, &samples_b);
    geomlib_bvh_free(&bvh_a);
    geomlib_bvh_free(&bvh_b);

    out->lower_bound = lower_a > lower_b ? lower_a : lower_b;
    out->upper_bound = upper_a > upper_b ? upper_a : upper_b;
    out->samples = samples_a + samples_b;
    return GEOMLIB_MESH_SUCCESS;
}
