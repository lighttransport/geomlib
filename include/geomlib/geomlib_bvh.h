/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * BVH (Bounding Volume Hierarchy) and Ray Tracing
 * Core functions implemented in C for maximum portability
 *
 * References:
 * - "Physically Based Rendering" pbr-book.org - BVH construction and SAH
 * - "A Survey on Bounding Volume Hierarchies for Ray Tracing" (2024)
 * - Jacco Bikker's BVH construction tutorials
 */

#ifndef GEOMLIB_BVH_H
#define GEOMLIB_BVH_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include "geomlib_linalg.h"
#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Ray Data Structure
 * ============================================================================ */

/**
 * Ray for ray tracing and ray casting queries
 */
typedef struct geomlib_ray {
    geomlib_vec3 origin;      /* Ray origin point */
    geomlib_vec3 direction;   /* Ray direction (should be normalized) */
    float t_min;              /* Minimum t parameter (default 0.0) */
    float t_max;              /* Maximum t parameter (default FLT_MAX) */
} geomlib_ray;

/**
 * Ray-surface intersection result
 */
typedef struct geomlib_ray_hit {
    bool hit;                 /* True if ray hit something */
    float t;                  /* Distance along ray to hit point */
    geomlib_vec3 position;    /* Hit position in world space */
    geomlib_vec3 normal;      /* Surface normal at hit point */
    geomlib_vec3 barycentric; /* Barycentric coordinates (for triangles) */
    uint32_t primitive_id;    /* Index of hit primitive (triangle) */
} geomlib_ray_hit;

/* ============================================================================
 * BVH Node Structure
 * ============================================================================ */

/**
 * BVH node (internal representation)
 *
 * Binary tree structure where each node contains an AABB.
 * Leaf nodes contain primitive indices, interior nodes have two children.
 */
typedef struct geomlib_bvh_node {
    geomlib_aabb bounds;      /* Bounding box for this node */

    /* Interior node: left and right child indices */
    /* Leaf node: primitive_offset and primitive_count */
    union {
        struct {
            uint32_t left_child;   /* Index of left child node */
            uint32_t right_child;  /* Index of right child node */
        } interior;
        struct {
            uint32_t primitive_offset;  /* Start index in primitive array */
            uint32_t primitive_count;   /* Number of primitives in leaf */
        } leaf;
    } data;

    uint8_t axis;             /* Split axis (0=x, 1=y, 2=z) */
    bool is_leaf;             /* True if this is a leaf node */
} geomlib_bvh_node;

/**
 * BVH structure
 *
 * Contains the complete hierarchy and reference to mesh data.
 */
typedef struct geomlib_bvh {
    geomlib_bvh_node *nodes;        /* Array of BVH nodes */
    uint32_t *primitive_indices;    /* Reordered primitive indices */
    size_t node_count;              /* Number of nodes */
    size_t node_capacity;           /* Allocated node capacity */
    size_t primitive_count;         /* Number of primitives */
    const geomlib_trimesh *mesh;    /* Reference to mesh (not owned) */
} geomlib_bvh;

/**
 * BVH build quality settings
 */
typedef enum geomlib_bvh_build_quality {
    GEOMLIB_BVH_QUALITY_LOW,      /* Fast build, lower quality */
    GEOMLIB_BVH_QUALITY_MEDIUM,   /* Balanced build quality */
    GEOMLIB_BVH_QUALITY_HIGH,     /* SAH with more buckets */
    GEOMLIB_BVH_QUALITY_FAST      /* Median split, fastest build */
} geomlib_bvh_build_quality;

/* ============================================================================
 * Ray-Primitive Intersection
 * ============================================================================ */

/**
 * Ray-triangle intersection using Möller-Trumbore algorithm
 *
 * Fast, watertight intersection test without computing plane equation.
 *
 * @param ray       Ray to test
 * @param v0        First triangle vertex
 * @param v1        Second triangle vertex
 * @param v2        Third triangle vertex
 * @param t_out     Output: distance along ray (if hit)
 * @param u_out     Output: barycentric u coordinate (if hit)
 * @param v_out     Output: barycentric v coordinate (if hit)
 * @return          True if ray intersects triangle
 */
bool geomlib_ray_triangle_intersect(
    const geomlib_ray *ray,
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    float *t_out,
    float *u_out,
    float *v_out
);

/**
 * Ray-AABB intersection test
 *
 * Optimized slab method for ray-box intersection.
 *
 * @param ray       Ray to test
 * @param aabb      Axis-aligned bounding box
 * @return          True if ray intersects AABB
 */
bool geomlib_ray_aabb_intersect(
    const geomlib_ray *ray,
    const geomlib_aabb *aabb
);

/* ============================================================================
 * BVH Construction
 * ============================================================================ */

/**
 * Build BVH from triangle mesh
 *
 * Constructs a binary bounding volume hierarchy using the specified
 * quality setting. Higher quality uses SAH (Surface Area Heuristic)
 * for better ray tracing performance but slower build times.
 *
 * @param bvh       BVH structure to initialize
 * @param mesh      Triangle mesh to build BVH for
 * @param quality   Build quality setting
 * @return          GEOMLIB_MESH_SUCCESS on success, error code otherwise
 */
int geomlib_bvh_build(
    geomlib_bvh *bvh,
    const geomlib_trimesh *mesh,
    geomlib_bvh_build_quality quality
);

/**
 * Free BVH resources
 *
 * @param bvh       BVH to free
 */
void geomlib_bvh_free(geomlib_bvh *bvh);

/**
 * Get BVH statistics
 *
 * @param bvh               BVH to analyze
 * @param max_depth_out     Output: maximum tree depth
 * @param leaf_count_out    Output: number of leaf nodes
 * @param avg_leaf_size_out Output: average primitives per leaf
 */
void geomlib_bvh_get_stats(
    const geomlib_bvh *bvh,
    uint32_t *max_depth_out,
    uint32_t *leaf_count_out,
    float *avg_leaf_size_out
);

/* ============================================================================
 * BVH Ray Queries
 * ============================================================================ */

/**
 * Trace ray through BVH (closest hit)
 *
 * Finds the closest intersection of ray with mesh geometry.
 *
 * @param bvh       BVH to query
 * @param ray       Ray to trace
 * @param hit       Output: hit information (if hit)
 * @return          True if ray hit anything
 */
bool geomlib_bvh_trace_ray(
    const geomlib_bvh *bvh,
    const geomlib_ray *ray,
    geomlib_ray_hit *hit
);

/**
 * Test ray occlusion (any hit)
 *
 * Faster than full trace when you only need to know if ray hits anything.
 * Returns immediately on first intersection.
 *
 * @param bvh       BVH to query
 * @param ray       Ray to test
 * @return          True if ray hits anything
 */
bool geomlib_bvh_test_occlusion(
    const geomlib_bvh *bvh,
    const geomlib_ray *ray
);

/* ============================================================================
 * BVH Point Queries
 * ============================================================================ */

/**
 * Find closest point on mesh to query point
 *
 * Uses BVH for acceleration. Much faster than brute-force for large meshes.
 *
 * @param bvh           BVH to query
 * @param query         Query point
 * @param closest_out   Output: closest point on mesh surface
 * @param distance_out  Output: distance to closest point (can be NULL)
 * @param primitive_id_out Output: triangle index of closest point (can be NULL)
 * @return              True if successful
 */
bool geomlib_bvh_closest_point(
    const geomlib_bvh *bvh,
    geomlib_vec3 query,
    geomlib_vec3 *closest_out,
    float *distance_out,
    uint32_t *primitive_id_out
);

/**
 * Find k nearest neighbors on mesh surface
 *
 * Returns up to k closest points on mesh surface to query point.
 * Results are sorted by distance (closest first).
 *
 * @param bvh           BVH to query
 * @param query         Query point
 * @param k             Number of neighbors to find
 * @param points_out    Output: array of k closest points
 * @param distances_out Output: array of k distances (can be NULL)
 * @param primitive_ids_out Output: array of k triangle indices (can be NULL)
 * @return              Number of neighbors found (may be less than k)
 */
size_t geomlib_bvh_knn_query(
    const geomlib_bvh *bvh,
    geomlib_vec3 query,
    size_t k,
    geomlib_vec3 *points_out,
    float *distances_out,
    uint32_t *primitive_ids_out
);

/**
 * Find all points within radius of query point
 *
 * Returns all mesh surface points within specified radius.
 *
 * @param bvh           BVH to query
 * @param query         Query point
 * @param radius        Search radius
 * @param max_results   Maximum number of results to return
 * @param points_out    Output: array of points within radius
 * @param distances_out Output: array of distances (can be NULL)
 * @param primitive_ids_out Output: array of triangle indices (can be NULL)
 * @return              Number of points found
 */
size_t geomlib_bvh_radius_query(
    const geomlib_bvh *bvh,
    geomlib_vec3 query,
    float radius,
    size_t max_results,
    geomlib_vec3 *points_out,
    float *distances_out,
    uint32_t *primitive_ids_out
);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * Create ray from origin and direction
 *
 * @param origin    Ray origin
 * @param direction Ray direction (will be normalized)
 * @return          Ray structure
 */
geomlib_ray geomlib_ray_create(geomlib_vec3 origin, geomlib_vec3 direction);

/**
 * Evaluate point along ray at parameter t
 *
 * @param ray       Ray
 * @param t         Parameter along ray
 * @return          Point at ray.origin + t * ray.direction
 */
geomlib_vec3 geomlib_ray_at(const geomlib_ray *ray, float t);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_BVH_H */
