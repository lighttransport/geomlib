/*
 * geomlib - Dihedral-angle region growing segmentation
 *
 * Grows connected face regions across edges whose dihedral angle is below a
 * threshold, producing a label per face.
 */

#ifndef GEOMLIB_REGION_GROWING_H
#define GEOMLIB_REGION_GROWING_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Label each triangle with a region id. Returns the number of regions.
 * `labels` must hold mesh->triangle_count entries. */
int geomlib_mesh_region_growing(const geomlib_trimesh *mesh,
                                float dihedral_threshold_deg,
                                uint32_t *labels, size_t *region_count);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_REGION_GROWING_H */
