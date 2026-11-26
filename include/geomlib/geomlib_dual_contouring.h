/*
 * geomlib - Dual contouring isosurface extraction
 *
 * Places one vertex per crossing grid cell by minimizing a quadric error
 * metric over the crossing edges' intersection planes, then connects the
 * cell vertices with quads around each crossing grid edge.
 */

#ifndef GEOMLIB_DUAL_CONTOURING_H
#define GEOMLIB_DUAL_CONTOURING_H

#include "geomlib_mesh.h"
#include "geomlib_sdf.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Extract the zero isosurface of a dense SDF with dual contouring. */
int geomlib_sdf_dual_contouring(const geomlib_sdf *sdf, geomlib_trimesh *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_DUAL_CONTOURING_H */
