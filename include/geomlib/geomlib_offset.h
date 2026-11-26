/*
 * geomlib - SDF-based triangle-mesh offsetting
 *
 * Evaluates a signed distance field around the input mesh and extracts the
 * level set at the requested offset: positive offsets grow the solid,
 * negative offsets shrink it.
 */

#ifndef GEOMLIB_OFFSET_H
#define GEOMLIB_OFFSET_H

#include "geomlib_sdf.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t resolution;   /* Grid resolution on the longest axis (default 48) */
    float padding;         /* Extra grid margin in voxels (default 6) */
} geomlib_offset_params;

geomlib_offset_params geomlib_offset_params_default(void);

/* Build an offset surface. Requires a closed, consistently oriented mesh. */
int geomlib_offset_mesh(const geomlib_trimesh *mesh, float offset,
                        const geomlib_offset_params *params,
                        geomlib_trimesh *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_OFFSET_H */
