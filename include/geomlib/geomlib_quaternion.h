/*
 * geomlib - Quaternion operations
 *
 * Unit-quaternion representation (w, x, y, z) with multiplication,
 * rotation-matrix conversion, normalization and spherical interpolation.
 */

#ifndef GEOMLIB_QUATERNION_H
#define GEOMLIB_QUATERNION_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double w, x, y, z; } geomlib_quat;

geomlib_quat geomlib_quat_identity(void);
geomlib_quat geomlib_quat_from_axis_angle(geomlib_vec3 axis, double angle);
geomlib_quat geomlib_quat_mul(geomlib_quat a, geomlib_quat b);
geomlib_quat geomlib_quat_conjugate(geomlib_quat q);
geomlib_quat geomlib_quat_normalize(geomlib_quat q);
double geomlib_quat_dot(geomlib_quat a, geomlib_quat b);
geomlib_mat4 geomlib_quat_to_mat4(geomlib_quat q);
geomlib_quat geomlib_quat_slerp(geomlib_quat a, geomlib_quat b, double t);
geomlib_vec3 geomlib_quat_rotate(geomlib_quat q, geomlib_vec3 v);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_QUATERNION_H */
