/*
 * geomlib - Marching tetrahedra isosurface extraction
 *
 * Splits each cube into six tetrahedra around a body diagonal and extracts
 * triangles from each tetrahedron's 16 sign configurations. Unlike marching
 * cubes there are no ambiguous faces. The output vertices are welded.
 */

#ifndef GEOMLIB_MARCHING_TETS_H
#define GEOMLIB_MARCHING_TETS_H

#include "geomlib_mesh.h"
#include "geomlib_sdf.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Extract the zero isosurface of a dense SDF using marching tetrahedra. */
int geomlib_sdf_marching_tets(const geomlib_sdf *sdf, geomlib_trimesh *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MARCHING_TETS_H */
