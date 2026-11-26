/*
 * geomlib - SDF-based mesh offsetting (see geomlib_offset.h)
 *
 * Clean-room implementation: sample signed distances on a padded grid around
 * the mesh, then extract the iso level equal to the requested offset.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_offset.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_offset_params geomlib_offset_params_default(void)
{
    geomlib_offset_params p;
    p.resolution = 48;
    p.padding = 6.0f;
    return p;
}

int geomlib_offset_mesh(const geomlib_trimesh *mesh, float offset,
                        const geomlib_offset_params *params,
                        geomlib_trimesh *out)
{
    geomlib_offset_params p;
    geomlib_sdf_params sp;
    geomlib_sdf field;
    geomlib_aabb box;
    geomlib_vec3 size;
    float extent;
    int rc;

    if (!mesh || !out || !mesh->vertices || !mesh->indices ||
        !mesh->triangle_count || !isfinite(offset)) return GEOMLIB_MESH_ERROR_INVALID;
    p = params ? *params : geomlib_offset_params_default();
    if (p.resolution < 16) p.resolution = 16;
    if (p.resolution > 128) p.resolution = 128;
    if (!(p.padding >= 2.0f)) p.padding = 2.0f;
    if (p.padding > 24.0f) p.padding = 24.0f;

    box = geomlib_trimesh_compute_aabb(mesh);
    size = geomlib_aabb_size(box);
    extent = fmaxf(size.x, fmaxf(size.y, size.z));
    if (!(extent > 1e-8f)) return GEOMLIB_MESH_ERROR_INVALID;

    sp = geomlib_sdf_params_default();
    sp.voxel_size = extent/(float)(p.resolution - 2.0f*p.padding);
    if (!(sp.voxel_size > 0.0f)) return GEOMLIB_MESH_ERROR_INVALID;
    sp.compute_sign = true;
    sp.use_bvh = true;
    sp.narrow_band = false;
    sp.sparse = false;
    sp.padding = p.padding;

    field = geomlib_mesh_to_sdf(mesh, &sp);
    if (!field.values) return GEOMLIB_MESH_ERROR_MEMORY;
    rc = geomlib_sdf_to_mesh(&field, offset, out);
    geomlib_sdf_free(&field);
    return rc;
}
