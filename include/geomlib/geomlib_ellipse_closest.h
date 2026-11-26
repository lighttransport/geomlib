/*
 * geomlib - Closest point on an ellipse
 *
 * Multi-start Newton solve of the stationary-distance condition.
 */

#ifndef GEOMLIB_ELLIPSE_CLOSEST_H
#define GEOMLIB_ELLIPSE_CLOSEST_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec2 point;   /* closest boundary point */
    float distance;       /* Euclidean distance to the query */
    float t;              /* ellipse parameter of the closest point */
} geomlib_ellipse_closest;

int geomlib_ellipse_closest_point(geomlib_vec2 q, geomlib_vec2 center,
                                  float a, float b, float angle,
                                  geomlib_ellipse_closest *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ELLIPSE_CLOSEST_H */
