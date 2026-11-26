/*
 * geomlib - Least-squares cylinder fitting
 *
 * Fits an infinite cylinder (axis direction, point on axis, radius) to 3D
 * points by minimizing the radial distance to the axis.
 */

#ifndef GEOMLIB_CYLINDER_FIT_H
#define GEOMLIB_CYLINDER_FIT_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 axis;      /* unit direction */
    geomlib_vec3 point;     /* a point on the axis */
    float radius;
    float rms;
} geomlib_cylinder_fit;

/* Fit a cylinder using `initial_axis` as a starting direction (need not be
 * unit). Returns GEOMLIB_MESH_SUCCESS on success. */
int geomlib_fit_cylinder(const geomlib_vec3 *points, size_t count,
                         geomlib_vec3 initial_axis,
                         geomlib_cylinder_fit *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_CYLINDER_FIT_H */
