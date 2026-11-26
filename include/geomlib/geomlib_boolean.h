/*
 * geomlib - Boolean operations on triangle meshes
 *
 * Surface-based booleans: both closed oriented surfaces are split along their
 * intersection curves, each fragment is classified against the other solid
 * with the generalized winding number, and the selected fragments are welded
 * into a consistently oriented output mesh.
 */

#ifndef GEOMLIB_BOOLEAN_H
#define GEOMLIB_BOOLEAN_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    GEOMLIB_BOOLEAN_UNION = 0,
    GEOMLIB_BOOLEAN_INTERSECTION,
    GEOMLIB_BOOLEAN_DIFFERENCE   /* A - B */
} geomlib_boolean_op;

/*
 * Compute the boolean operation on two closed, consistently oriented meshes.
 * `weld_eps` is the vertex-merging tolerance (0 selects 1e-5 * bbox diagonal).
 * The output is a newly initialized mesh owned by the caller.
 */
int geomlib_mesh_boolean(const geomlib_trimesh *a,
                         const geomlib_trimesh *b,
                         geomlib_boolean_op op,
                         float weld_eps,
                         geomlib_trimesh *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_BOOLEAN_H */
