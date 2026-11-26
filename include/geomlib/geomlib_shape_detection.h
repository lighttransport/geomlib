/*
 * geomlib - RANSAC shape detection for point clouds
 *
 * Repeatedly hypotheses geometric primitives (planes, spheres, cylinders)
 * from minimal point samples, scores them by inlier count, and removes the
 * best primitive's inliers until the requested number of shapes is found.
 */

#ifndef GEOMLIB_SHAPE_DETECTION_H
#define GEOMLIB_SHAPE_DETECTION_H

#include "geomlib_surface.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    GEOMLIB_SHAPE_PLANE = 0,
    GEOMLIB_SHAPE_SPHERE,
    GEOMLIB_SHAPE_CYLINDER
} geomlib_shape_type;

typedef struct {
    geomlib_shape_type type;
    float threshold;        /* Distance threshold for inliers */
    geomlib_vec3 center;    /* Sphere center / cylinder axis origin / plane point */
    geomlib_vec3 normal;    /* Plane normal / cylinder axis direction */
    float radius;           /* Sphere/cylinder radius */
    size_t inlier_count;
    uint32_t *inliers;      /* Indices into the input cloud */
} geomlib_detected_shape;

geomlib_detected_shape geomlib_detected_shape_create(void);
void geomlib_detected_shape_free(geomlib_detected_shape *s);

/*
 * Detect up to `max_shapes` primitives in `cloud`. Returns the number found
 * in `out_count`; the caller owns the returned array and frees each shape.
 */
geomlib_detected_shape *geomlib_detect_shapes_ransac(
    const geomlib_pointcloud *cloud,
    geomlib_shape_type type,
    float threshold,
    size_t max_shapes,
    uint32_t seed,
    size_t *out_count);

void geomlib_detected_shapes_free(geomlib_detected_shape *shapes, size_t count);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SHAPE_DETECTION_H */
