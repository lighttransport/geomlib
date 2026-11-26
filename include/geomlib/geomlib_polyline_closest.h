/*
 * geomlib - Closest point on a polyline
 *
 * Projects a query onto every segment and keeps the nearest result.
 */

#ifndef GEOMLIB_POLYLINE_CLOSEST_H
#define GEOMLIB_POLYLINE_CLOSEST_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 point;   /* closest point on the polyline */
    float distance;       /* Euclidean distance to the query */
    size_t segment;       /* index of the closest segment */
    float t;              /* parameter in [0,1] along that segment */
} geomlib_polyline_closest;

int geomlib_polyline_closest_point(const geomlib_vec3 *points, size_t n,
                                   geomlib_vec3 q,
                                   geomlib_polyline_closest *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_POLYLINE_CLOSEST_H */
