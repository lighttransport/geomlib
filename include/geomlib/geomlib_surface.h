/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Surface Reconstruction from Point Clouds
 * Fast algorithms optimized for realtime use and Gaussian splats
 *
 * References:
 * - "A Simple, Fast, and Effective Polygon Reduction Algorithm" (1998)
 * - "The Ball-Pivoting Algorithm for Surface Reconstruction" (1999)
 * - "Poisson Surface Reconstruction" (2006)
 * - "Fast and Robust Surface Reconstruction from Range Data" (2002)
 */

#ifndef GEOMLIB_SURFACE_H
#define GEOMLIB_SURFACE_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include "geomlib_linalg.h"
#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Point Cloud Data Structure
 * ============================================================================ */

/**
 * Point cloud structure
 *
 * Stores unorganized 3D points with optional normals and attributes.
 * Used as input for surface reconstruction algorithms.
 */
typedef struct geomlib_pointcloud {
    geomlib_vec3 *points;        /* Point positions */
    geomlib_vec3 *normals;       /* Point normals (optional, can be NULL) */
    float *confidences;          /* Point confidence/weight (optional, can be NULL) */
    size_t point_count;          /* Number of points */
    size_t point_capacity;       /* Allocated capacity */
} geomlib_pointcloud;

/* ============================================================================
 * Point Cloud Management
 * ============================================================================ */

/**
 * Initialize point cloud
 *
 * @param pc                Point cloud to initialize
 * @param initial_capacity  Initial capacity for points
 * @return                  GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_pointcloud_init(geomlib_pointcloud *pc, size_t initial_capacity);

/**
 * Free point cloud resources
 *
 * @param pc    Point cloud to free
 */
void geomlib_pointcloud_free(geomlib_pointcloud *pc);

/**
 * Add point to point cloud
 *
 * @param pc            Point cloud
 * @param position      Point position
 * @param normal        Point normal (can be NULL)
 * @param confidence    Point confidence/weight (0-1, use 1.0 for uniform)
 * @return              Index of added point
 */
size_t geomlib_pointcloud_add_point(
    geomlib_pointcloud *pc,
    geomlib_vec3 position,
    geomlib_vec3 *normal,
    float confidence
);

/**
 * Estimate normals from k-nearest neighbors
 *
 * Uses PCA to estimate normal direction from local neighborhood.
 * Orients normals consistently using minimum spanning tree.
 *
 * @param pc    Point cloud
 * @param k     Number of neighbors to use (typically 10-30)
 * @return      GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_pointcloud_estimate_normals(geomlib_pointcloud *pc, int k);

/**
 * Downsample point cloud using voxel grid
 *
 * Reduces point count by averaging points within voxel cells.
 *
 * @param pc            Input point cloud
 * @param voxel_size    Size of voxel cells
 * @param output        Output downsampled point cloud
 * @return              GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_pointcloud_downsample(
    const geomlib_pointcloud *pc,
    float voxel_size,
    geomlib_pointcloud *output
);

/* ============================================================================
 * Surface Reconstruction Algorithms
 * ============================================================================ */

/**
 * Surface reconstruction method
 */
typedef enum geomlib_surface_method {
    GEOMLIB_SURFACE_GREEDY_PROJECTION,  /* Fast greedy triangulation */
    GEOMLIB_SURFACE_BALL_PIVOTING,      /* Ball pivoting algorithm */
    GEOMLIB_SURFACE_MARCHING_CUBES,     /* Voxel-based marching cubes */
    GEOMLIB_SURFACE_ALPHA_SHAPE         /* Alpha shape reconstruction */
} geomlib_surface_method;

/**
 * Surface reconstruction parameters
 */
typedef struct geomlib_surface_params {
    geomlib_surface_method method;  /* Reconstruction method */

    /* Common parameters */
    float search_radius;            /* Neighborhood search radius */
    int max_neighbors;              /* Maximum neighbors to consider */

    /* Greedy projection parameters */
    float max_edge_length;          /* Maximum allowed edge length */
    float max_surface_angle;        /* Maximum angle between surface normals (degrees) */

    /* Ball pivoting parameters */
    float ball_radius;              /* Radius of ball for pivoting */

    /* Marching cubes parameters */
    float voxel_size;               /* Size of voxel grid cells */
    float iso_value;                /* Iso-surface value (typically 0.0) */

    /* Alpha shape parameters */
    float alpha;                    /* Alpha parameter for alpha shapes */
} geomlib_surface_params;

/**
 * Create default surface reconstruction parameters
 *
 * @param method    Reconstruction method
 * @return          Default parameters for method
 */
geomlib_surface_params geomlib_surface_params_default(geomlib_surface_method method);

/**
 * Reconstruct surface from point cloud
 *
 * Fast surface reconstruction optimized for realtime use.
 * Produces low-poly proxy mesh suitable for collision detection and preview.
 *
 * @param pc        Input point cloud (must have normals for some methods)
 * @param params    Reconstruction parameters
 * @param output    Output triangle mesh
 * @return          GEOMLIB_MESH_SUCCESS on success, error code otherwise
 */
int geomlib_reconstruct_surface(
    const geomlib_pointcloud *pc,
    const geomlib_surface_params *params,
    geomlib_trimesh *output
);

/**
 * Quick surface reconstruction with automatic parameters
 *
 * Analyzes point cloud and chooses optimal parameters automatically.
 * Best for quick previews and prototyping.
 *
 * @param pc            Input point cloud
 * @param target_faces  Target number of faces (0 = automatic)
 * @param output        Output triangle mesh
 * @return              GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_reconstruct_surface_auto(
    const geomlib_pointcloud *pc,
    size_t target_faces,
    geomlib_trimesh *output
);

/* ============================================================================
 * Gaussian Splat Conversion
 * ============================================================================ */

/**
 * Convert Gaussian splats to point cloud
 *
 * Extracts point positions from Gaussian splat data.
 * Centers and scales are used to generate representative points.
 *
 * @param positions     Gaussian centers (xyz)
 * @param count         Number of Gaussians
 * @param output        Output point cloud
 * @return              GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_gaussians_to_pointcloud(
    const float *positions,
    size_t count,
    geomlib_pointcloud *output
);

/**
 * Sample surface points from Gaussian splats
 *
 * Generates additional points by sampling from Gaussian distributions.
 * Useful for denser reconstruction.
 *
 * @param positions     Gaussian centers (xyz)
 * @param scales        Gaussian scales (3 per Gaussian, can be NULL for uniform)
 * @param rotations     Gaussian rotations (4 quaternion components, can be NULL)
 * @param count         Number of Gaussians
 * @param samples_per   Number of samples per Gaussian
 * @param output        Output point cloud
 * @return              GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_gaussians_sample_points(
    const float *positions,
    const float *scales,
    const float *rotations,
    size_t count,
    int samples_per,
    geomlib_pointcloud *output
);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * Compute point cloud bounding box
 *
 * @param pc        Point cloud
 * @return          Bounding box
 */
geomlib_aabb geomlib_pointcloud_bounds(const geomlib_pointcloud *pc);

/**
 * Get point cloud statistics
 *
 * @param pc                Point cloud
 * @param avg_spacing_out   Output: average point spacing
 * @param density_out       Output: point density (points per unit volume)
 */
void geomlib_pointcloud_get_stats(
    const geomlib_pointcloud *pc,
    float *avg_spacing_out,
    float *density_out
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SURFACE_H */
