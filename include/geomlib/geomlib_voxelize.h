/*
 * geomlib - Triangle-mesh voxelization
 *
 * Marks voxels intersected by the surface (conservative SAT overlap) and,
 * for closed meshes, voxels whose centres lie inside the surface (BVH ray
 * parity).
 */

#ifndef GEOMLIB_VOXELIZE_H
#define GEOMLIB_VOXELIZE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t nx, ny, nz;
    geomlib_vec3 origin;
    float voxel_size;
    uint8_t *surface;         /* nx*ny*nz, 1 when the voxel meets the surface */
    uint8_t *volume;          /* nx*ny*nz, 1 when the centre is inside */
    size_t surface_count;
    size_t volume_count;
} geomlib_voxel_grid;

void geomlib_voxel_grid_free(geomlib_voxel_grid *g);

/* Rasterize `mesh` into an nx*ny*nz grid spanning the mesh bounds with a
 * 5% margin. The returned grid owns both label arrays (empty on failure). */
geomlib_voxel_grid geomlib_voxelize_mesh(const geomlib_trimesh *mesh,
                                         uint32_t nx, uint32_t ny, uint32_t nz);

/* Linear index of a voxel, or SIZE_MAX when out of range. */
size_t geomlib_voxel_grid_index(const geomlib_voxel_grid *g,
                                uint32_t x, uint32_t y, uint32_t z);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_VOXELIZE_H */
