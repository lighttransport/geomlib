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

/* Opaque sparse narrow-band storage. NULL for dense grids. */
struct geomlib_sdf_band;

/* Signed distance field */
typedef struct geomlib_sdf {
    float *values;                  /* Distance values (linearized 3D array), NULL when sparse */

    geomlib_vec3 origin;            /* Grid origin (min corner) */
    geomlib_vec3 extent;            /* Grid extent (size) */
    float voxel_size;               /* Size of each voxel */

    uint32_t resolution[3];         /* Grid resolution (nx, ny, nz) */
    size_t total_voxels;            /* nx * ny * nz */

    bool is_narrow_band;            /* Is this a narrow-band SDF? */
    float narrow_band_width;        /* Band width in voxels (0 = full SDF) */

    bool is_sparse;                 /* True when only defined samples are stored */
    size_t stored_voxels;           /* Number of defined samples */

    geomlib_aabb bounds;            /* Bounding box */

    struct geomlib_sdf_band *band;  /* Sparse storage when is_sparse */
} geomlib_sdf;

/* Distance evaluation strategy */
typedef enum geomlib_sdf_algorithm {
    GEOMLIB_SDF_ALGORITHM_BRUTE_FORCE = 0,   /* Exact per-voxel distance (default) */
    GEOMLIB_SDF_ALGORITHM_FAST_SWEEPING = 1  /* Exact band seed, then eikonal sweeping */
} geomlib_sdf_algorithm;

/* SDF generation parameters */
typedef struct geomlib_sdf_params {
    float voxel_size;               /* Size of each voxel (0 = auto) */
    uint32_t resolution;            /* Target resolution (used if voxel_size = 0) */

    bool narrow_band;               /* Generate narrow-band SDF */
    float band_width;               /* Band width in voxels (default: 3.0) */

    bool compute_sign;              /* Compute inside/outside sign */
    bool use_bvh;                   /* Use BVH for acceleration (recommended) */

    float padding;                  /* Extra padding around mesh (in voxels) */

    bool sparse;                    /* Store only the band; requires narrow_band */
    geomlib_sdf_algorithm algorithm; /* Distance evaluation strategy */
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

/* True when the sample at (x, y, z) has a stored value. Out-of-range is false. */
bool geomlib_sdf_is_defined(
    const geomlib_sdf *sdf,
    uint32_t x,
    uint32_t y,
    uint32_t z
);

/* Number of stored samples. Equals total_voxels for dense grids. */
size_t geomlib_sdf_stored_voxels(const geomlib_sdf *sdf);

/* Convert a sparse band into a dense grid. Undefined samples are clamped to the
 * band limit with the sign of the nearest defined neighbour. Returns
 * GEOMLIB_MESH_SUCCESS on an already-dense grid. Ownership stays with the caller. */
int geomlib_sdf_make_dense(geomlib_sdf *sdf);

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

/* Query by trilinear interpolation between voxel centers. Returns FLT_MAX outside that domain. */
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

/* Query gradient using central differences, with one-sided differences at the boundary. */
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

/* Extract an indexed isosurface from voxel-center samples. On success output owns a new mesh. */
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

/* Extract the zero-level set and place vertices on reference sharp edges.
 * Vertices within 1.5 voxels of a reference feature segment snap onto it, which
 * preserves creases. When snap_to_surface is true other vertices project onto the
 * closest reference surface point. reference and its features may be NULL, in
 * which case this matches geomlib_sdf_extract_surface. */
int geomlib_sdf_extract_surface_features(
    const geomlib_sdf *sdf,
    const geomlib_trimesh *reference,
    float sharp_angle_deg,
    bool snap_to_surface,
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

/* Resample the longest axis to new_resolution, preserving origin and uniform spacing.
 * Dense input only; materialize sparse grids with geomlib_sdf_make_dense first.
 * The returned grid owns new storage. An empty result indicates failure. */
geomlib_sdf geomlib_sdf_resize(
    const geomlib_sdf *sdf,
    uint32_t new_resolution
);

/* Copy samples whose centers lie inside region, preserving voxel spacing.
 * The returned grid owns new storage. An empty result indicates no samples or failure. */
geomlib_sdf geomlib_sdf_extract_region(
    const geomlib_sdf *sdf,
    geomlib_aabb region
);

/* ============================================================================
 * File I/O
 * ============================================================================ */

/* Save SDF to file (binary format) */
int geomlib_sdf_save(
    const geomlib_sdf *sdf,
    const char *filename
);

/* Load an SDF1 grid into an empty output. The output is assigned only on success. */
int geomlib_sdf_load(
    geomlib_sdf *sdf,
    const char *filename
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SDF_H */
