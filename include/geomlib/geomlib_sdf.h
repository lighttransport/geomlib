/*
 * geomlib - Signed Distance Field (SDF) Generation
 * Convert meshes to volumetric distance fields and query distance/gradient
 */

#ifndef GEOMLIB_SDF_H
#define GEOMLIB_SDF_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include "geomlib_bvh.h"
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * SDF Data Structures
 * ============================================================================ */

/* Signed distance field */
typedef struct geomlib_sdf {
    float *values;                  /* Distance values (linearized 3D array) */

    geomlib_vec3 origin;            /* Grid origin (min corner) */
    geomlib_vec3 extent;            /* Grid extent (size) */
    float voxel_size;               /* Size of each voxel */

    uint32_t resolution[3];         /* Grid resolution (nx, ny, nz) */
    size_t total_voxels;            /* nx * ny * nz */

    bool is_narrow_band;            /* Is this a narrow-band SDF? */
    float narrow_band_width;        /* Band width in voxels (0 = full SDF) */

    geomlib_aabb bounds;            /* Bounding box */
} geomlib_sdf;

/* SDF generation parameters */
typedef struct geomlib_sdf_params {
    float voxel_size;               /* Size of each voxel (0 = auto) */
    uint32_t resolution;            /* Target resolution (used if voxel_size = 0) */

    bool narrow_band;               /* Generate narrow-band SDF */
    float band_width;               /* Band width in voxels (default: 3.0) */

    bool compute_sign;              /* Compute inside/outside sign */
    bool use_bvh;                   /* Use BVH for acceleration (recommended) */

    float padding;                  /* Extra padding around mesh (in voxels) */
} geomlib_sdf_params;

/* ============================================================================
 * SDF Generation
 * ============================================================================ */

/* Default SDF parameters */
geomlib_sdf_params geomlib_sdf_params_default(void);

/* Generate SDF from triangle mesh */
geomlib_sdf geomlib_mesh_to_sdf(
    const geomlib_trimesh *mesh,
    const geomlib_sdf_params *params
);

/* Generate SDF with explicit voxel size */
geomlib_sdf geomlib_mesh_to_sdf_voxel_size(
    const geomlib_trimesh *mesh,
    float voxel_size
);

/* Generate SDF with target resolution */
geomlib_sdf geomlib_mesh_to_sdf_resolution(
    const geomlib_trimesh *mesh,
    uint32_t resolution
);

/* Free SDF data */
void geomlib_sdf_free(geomlib_sdf *sdf);

/* ============================================================================
 * SDF Queries
 * ============================================================================ */

/* Get voxel index from 3D coordinates */
size_t geomlib_sdf_get_index(
    const geomlib_sdf *sdf,
    uint32_t x,
    uint32_t y,
    uint32_t z
);

/* Get voxel coordinates from world position */
bool geomlib_sdf_world_to_voxel(
    const geomlib_sdf *sdf,
    geomlib_vec3 world_pos,
    int32_t *x,
    int32_t *y,
    int32_t *z
);

/* Query distance at world position (trilinear interpolation) */
float geomlib_sdf_query_distance(
    const geomlib_sdf *sdf,
    geomlib_vec3 point
);

/* Query distance at voxel coordinates (no interpolation) */
float geomlib_sdf_query_distance_voxel(
    const geomlib_sdf *sdf,
    uint32_t x,
    uint32_t y,
    uint32_t z
);

/* Query gradient at world position (central differences) */
geomlib_vec3 geomlib_sdf_query_gradient(
    const geomlib_sdf *sdf,
    geomlib_vec3 point
);

/* Query normal at world position (normalized gradient) */
geomlib_vec3 geomlib_sdf_query_normal(
    const geomlib_sdf *sdf,
    geomlib_vec3 point
);

/* Check if point is inside mesh (distance < 0) */
bool geomlib_sdf_is_inside(
    const geomlib_sdf *sdf,
    geomlib_vec3 point
);

/* ============================================================================
 * SDF to Mesh Conversion (Marching Cubes)
 * ============================================================================ */

/* Extract isosurface using Marching Cubes */
int geomlib_sdf_to_mesh(
    const geomlib_sdf *sdf,
    float iso_value,
    geomlib_trimesh *output
);

/* Extract zero-level set (iso_value = 0) */
int geomlib_sdf_extract_surface(
    const geomlib_sdf *sdf,
    geomlib_trimesh *output
);

/* ============================================================================
 * Distance Field Computation
 * ============================================================================ */

/* Compute unsigned distance to mesh at point */
float geomlib_compute_unsigned_distance(
    const geomlib_trimesh *mesh,
    const geomlib_bvh *bvh,
    geomlib_vec3 point
);

/* Compute signed distance to mesh at point */
float geomlib_compute_signed_distance(
    const geomlib_trimesh *mesh,
    const geomlib_bvh *bvh,
    geomlib_vec3 point
);

/* Find closest point on mesh to query point */
geomlib_vec3 geomlib_find_closest_point_on_mesh(
    const geomlib_trimesh *mesh,
    const geomlib_bvh *bvh,
    geomlib_vec3 point,
    float *distance
);

/* ============================================================================
 * Sign Determination
 * ============================================================================ */

/* Determine if point is inside mesh (ray casting) */
bool geomlib_is_point_inside_mesh(
    const geomlib_trimesh *mesh,
    const geomlib_bvh *bvh,
    geomlib_vec3 point
);

/* Determine sign using generalized winding number */
float geomlib_compute_winding_number(
    const geomlib_trimesh *mesh,
    geomlib_vec3 point
);

/* ============================================================================
 * SDF Operations
 * ============================================================================ */

/* Smooth SDF with Gaussian filter */
int geomlib_sdf_smooth(
    geomlib_sdf *sdf,
    float sigma
);

/* Dilate SDF (expand surface) */
int geomlib_sdf_dilate(
    geomlib_sdf *sdf,
    float amount
);

/* Erode SDF (shrink surface) */
int geomlib_sdf_erode(
    geomlib_sdf *sdf,
    float amount
);

/* Blend two SDFs */
geomlib_sdf geomlib_sdf_blend(
    const geomlib_sdf *sdf1,
    const geomlib_sdf *sdf2,
    float blend_factor
);

/* ============================================================================
 * SDF Utilities
 * ============================================================================ */

/* Compute optimal voxel size for mesh */
float geomlib_compute_optimal_voxel_size(
    const geomlib_trimesh *mesh,
    uint32_t target_resolution
);

/* Compute SDF statistics */
typedef struct geomlib_sdf_stats {
    float min_distance;
    float max_distance;
    float avg_distance;

    size_t inside_voxels;
    size_t outside_voxels;
    size_t surface_voxels;      /* Voxels near surface */

    size_t memory_bytes;
} geomlib_sdf_stats;

geomlib_sdf_stats geomlib_compute_sdf_stats(const geomlib_sdf *sdf);

/* Resize SDF to new resolution */
geomlib_sdf geomlib_sdf_resize(
    const geomlib_sdf *sdf,
    uint32_t new_resolution
);

/* Copy SDF region */
geomlib_sdf geomlib_sdf_extract_region(
    const geomlib_sdf *sdf,
    geomlib_aabb region
);

/* ============================================================================
 * File I/O (Future)
 * ============================================================================ */

/* Save SDF to file (binary format) */
int geomlib_sdf_save(
    const geomlib_sdf *sdf,
    const char *filename
);

/* Load SDF from file */
int geomlib_sdf_load(
    geomlib_sdf *sdf,
    const char *filename
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SDF_H */
