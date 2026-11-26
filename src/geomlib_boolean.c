/*
 * geomlib - Boolean operations on triangle meshes (see geomlib_boolean.h)
 *
 * Grid-based SDF booleans. Both solids are sampled on one shared grid built
 * from the union of their bounds; min / max / max(dA,-dB) combine the signed
 * distances and the existing marching-cubes path extracts the result surface.
 * `voxel_hint` <= 0 selects max_extent / 64.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_boolean.h"
#include "geomlib/geomlib_sdf.h"
#include "geomlib/geomlib_bvh.h"

#include <math.h>
#include <stdlib.h>

int geomlib_mesh_boolean(const geomlib_trimesh *a,
                         const geomlib_trimesh *b,
                         geomlib_boolean_op op,
                         float voxel_hint,
                         geomlib_trimesh *out)
{
    geomlib_aabb ba, bb, both;
    geomlib_vec3 size, center;
    geomlib_sdf_params params;
    geomlib_sdf grid;
    geomlib_bvh bvh_a, bvh_b;
    float max_extent, h;
    uint32_t nx, ny, nz;
    size_t x, y, z;
    int rc;

    if (!a || !b || !out || !a->vertices || !b->vertices ||
        a->triangle_count == 0 || b->triangle_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;

    ba = geomlib_trimesh_compute_aabb(a);
    bb = geomlib_trimesh_compute_aabb(b);
    both = geomlib_aabb_merge(ba, bb);
    size = geomlib_aabb_size(both);
    center = geomlib_aabb_center(both);
    max_extent = fmaxf(size.x, fmaxf(size.y, size.z));
    h = (voxel_hint > 0.0f) ? voxel_hint : max_extent / 32.0f;
    if (!(h > 0.0f)) h = 1e-3f;
    nx = (uint32_t)(size.x / h) + 5;
    ny = (uint32_t)(size.y / h) + 5;
    nz = (uint32_t)(size.z / h) + 5;
    if (nx < 8) nx = 8;
    if (ny < 8) ny = 8;
    if (nz < 8) nz = 8;
    if (nx > 256 || ny > 256 || nz > 256) return GEOMLIB_MESH_ERROR_INVALID;

    params = geomlib_sdf_params_default();
    params.voxel_size = h;
    params.resolution = (nx > ny) ? (nx > nz ? nx : nz) : (ny > nz ? ny : nz);
    params.compute_sign = true;
    params.use_bvh = true;
    params.algorithm = GEOMLIB_SDF_ALGORITHM_BRUTE_FORCE;

    /* Build a one-triangle placeholder grid just to obtain a valid dense SDF
     * with the exact target origin/resolution/voxel size. */
    {
        geomlib_trimesh dummy = {0};
        geomlib_create_box(&dummy, size.x, size.y, size.z);
        for (size_t i = 0; i < dummy.vertex_count; ++i)
            dummy.vertices[i] = geomlib_vec3_add(dummy.vertices[i], center);
        params.voxel_size = h;
        /* create_box spans +/- size/2, so after translation the dummy box has
         * exactly the combined bounds: the generated grid is the shared grid. */
        grid = geomlib_mesh_to_sdf(&dummy, &params);
        geomlib_trimesh_free(&dummy);
    }
    if (!grid.values) return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_bvh_build(&bvh_a, a, GEOMLIB_BVH_QUALITY_MEDIUM);
    geomlib_bvh_build(&bvh_b, b, GEOMLIB_BVH_QUALITY_MEDIUM);

    for (z = 0; z < grid.resolution[2]; ++z) {
        for (y = 0; y < grid.resolution[1]; ++y) {
            for (x = 0; x < grid.resolution[0]; ++x) {
                geomlib_vec3 p = {
                    grid.origin.x + ((float)x + 0.5f) * h,
                    grid.origin.y + ((float)y + 0.5f) * h,
                    grid.origin.z + ((float)z + 0.5f) * h
                };
                float da = geomlib_compute_signed_distance(a, &bvh_a, p);
                float db = geomlib_compute_signed_distance(b, &bvh_b, p);
                size_t idx = geomlib_sdf_get_index(&grid, (uint32_t)x, (uint32_t)y, (uint32_t)z);
                switch (op) {
                    case GEOMLIB_BOOLEAN_UNION:        grid.values[idx] = fminf(da, db); break;
                    case GEOMLIB_BOOLEAN_INTERSECTION: grid.values[idx] = fmaxf(da, db); break;
                    case GEOMLIB_BOOLEAN_DIFFERENCE:   grid.values[idx] = fmaxf(da, -db); break;
                    default: break;
                }
            }
        }
    }

    rc = geomlib_sdf_extract_surface(&grid, out);
    geomlib_bvh_free(&bvh_a);
    geomlib_bvh_free(&bvh_b);
    geomlib_sdf_free(&grid);
    return rc;
}
