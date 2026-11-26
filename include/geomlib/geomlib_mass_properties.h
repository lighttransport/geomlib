/*
 * geomlib - Mass properties of a closed triangle mesh
 *
 * Computes volume, centre of mass and inertia tensor assuming unit density,
 * using signed tetrahedron decomposition against the origin.
 */

#ifndef GEOMLIB_MASS_PROPERTIES_H
#define GEOMLIB_MASS_PROPERTIES_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float volume;
    geomlib_vec3 center;
    geomlib_mat3 inertia;      /* column-major, about the centre of mass */
} geomlib_mass_properties;

/* Requires a closed, consistently oriented mesh with positive volume. */
int geomlib_compute_mass_properties(const geomlib_trimesh *mesh,
                                    geomlib_mass_properties *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MASS_PROPERTIES_H */
