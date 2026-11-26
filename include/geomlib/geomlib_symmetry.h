/*
 * geomlib - Reflection symmetry detection
 *
 * Tests the three principal planes through the centroid and reports the one
 * whose reflected point set best matches the original surface.
 */

#ifndef GEOMLIB_SYMMETRY_H
#define GEOMLIB_SYMMETRY_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    geomlib_vec3 plane_point;
    geomlib_vec3 plane_normal;
    float rms_error;
} geomlib_symmetry_result;

/* Detect the best reflection plane. A lower rms_error means more symmetric. */
int geomlib_detect_symmetry(const geomlib_trimesh *mesh,
                            geomlib_symmetry_result *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SYMMETRY_H */
