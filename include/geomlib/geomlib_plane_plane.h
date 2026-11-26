/*
 * geomlib - Plane / plane intersection
 *
 * Computes the line where two infinite planes meet.
 */

#ifndef GEOMLIB_PLANE_PLANE_H
#define GEOMLIB_PLANE_PLANE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 point;      /* a point on the intersection line */
    geomlib_vec3 direction;  /* unit direction of the line */
    int parallel;            /* 1 when the planes are parallel or coincident */
} geomlib_plane_line;

/* Returns 1 when the planes meet in a line, else 0 (parallel/coincident). */
int geomlib_plane_plane_intersection(geomlib_vec3 p0, geomlib_vec3 n0,
                                     geomlib_vec3 p1, geomlib_vec3 n1,
                                     geomlib_plane_line *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_PLANE_PLANE_H */
