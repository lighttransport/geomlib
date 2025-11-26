/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Collision Detection and Mesh Reduction
 * Convex hull generation, mesh simplification, proxy mesh generation
 *
 * References:
 * - "Real-Time Collision Detection" by Christer Ericson
 * - Quickhull algorithm (Barber, Dobkin, Huhdanpaa 1996)
 * - "Surface Simplification Using Quadric Error Metrics" (Garland & Heckbert 1997)
 */

#ifndef GEOMLIB_COLLISION_H
#define GEOMLIB_COLLISION_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include "geomlib_linalg.h"
#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Oriented Bounding Box (OBB)
 * ============================================================================ */

/**
 * Oriented bounding box
 * More accurate than AABB for arbitrarily oriented objects
 */
typedef struct geomlib_obb {
    geomlib_vec3 center;        /* Center of OBB */
    geomlib_vec3 axes[3];       /* Local axes (orthonormal) */
    geomlib_vec3 half_extents;  /* Half-sizes along each axis */
} geomlib_obb;

/**
 * Compute oriented bounding box for mesh using PCA
 *
 * Uses Principal Component Analysis to find optimal orientation
 * that minimizes bounding volume.
 *
 * @param mesh      Input mesh
 * @param obb_out   Output OBB
 * @return          GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_compute_obb(const geomlib_trimesh *mesh, geomlib_obb *obb_out);

/**
 * Compute OBB volume
 *
 * @param obb       OBB
 * @return          Volume
 */
float geomlib_obb_volume(const geomlib_obb *obb);

/**
 * Test if point is inside OBB
 *
 * @param obb       OBB
 * @param point     Point to test
 * @return          True if point is inside
 */
bool geomlib_obb_contains_point(const geomlib_obb *obb, geomlib_vec3 point);

/* ============================================================================
 * Convex Hull Generation
 * ============================================================================ */

/**
 * Compute convex hull of mesh using Quickhull algorithm
 *
 * Generates a convex hull mesh that encloses all vertices.
 * Result is suitable for collision detection proxies.
 *
 * @param mesh          Input mesh
 * @param hull_out      Output convex hull mesh (must be initialized)
 * @return              GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_compute_convex_hull(const geomlib_trimesh *mesh, geomlib_trimesh *hull_out);

/**
 * Test if mesh is convex
 *
 * @param mesh      Mesh to test
 * @param epsilon   Tolerance for convexity test
 * @return          True if mesh is convex
 */
bool geomlib_is_convex(const geomlib_trimesh *mesh, float epsilon);

/* ============================================================================
 * Mesh Simplification
 * ============================================================================ */

/**
 * Simplify mesh using edge collapse with quadric error metrics
 *
 * Reduces triangle count while preserving mesh shape.
 * Uses Garland-Heckbert quadric error metric to determine
 * which edges to collapse.
 *
 * @param mesh              Input mesh (will be modified in-place)
 * @param target_triangles  Target triangle count (must be > 4)
 * @param preserve_boundary Preserve boundary edges (1) or not (0)
 * @return                  GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_simplify_mesh(geomlib_trimesh *mesh, size_t target_triangles, int preserve_boundary);

/**
 * Simplify mesh by target reduction ratio
 *
 * @param mesh              Input mesh (will be modified in-place)
 * @param reduction_ratio   Reduction ratio (0.0 = no reduction, 1.0 = maximum reduction)
 * @param preserve_boundary Preserve boundary edges
 * @return                  GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_simplify_mesh_ratio(geomlib_trimesh *mesh, float reduction_ratio, int preserve_boundary);

/* ============================================================================
 * Collision Proxy Mesh Generation
 * ============================================================================ */

/**
 * Proxy mesh generation strategy
 */
typedef enum geomlib_proxy_strategy {
    GEOMLIB_PROXY_CONVEX_HULL,      /* Convex hull (fastest collision, lowest accuracy) */
    GEOMLIB_PROXY_SIMPLIFIED,       /* Simplified mesh (balanced) */
    GEOMLIB_PROXY_OBB,              /* Oriented bounding box (very fast, low accuracy) */
    GEOMLIB_PROXY_AABB              /* Axis-aligned bounding box (fastest, lowest accuracy) */
} geomlib_proxy_strategy;

/**
 * Generate collision proxy mesh
 *
 * Creates a simplified mesh suitable for collision detection.
 * The generated mesh has fewer triangles and is optimized for
 * fast collision queries.
 *
 * @param mesh              Input high-resolution mesh
 * @param strategy          Proxy generation strategy
 * @param target_triangles  Target triangle count (for SIMPLIFIED strategy)
 * @param proxy_out         Output proxy mesh (must be initialized)
 * @return                  GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_generate_collision_proxy(
    const geomlib_trimesh *mesh,
    geomlib_proxy_strategy strategy,
    size_t target_triangles,
    geomlib_trimesh *proxy_out
);

/**
 * Generate bounding box mesh (6 faces, 12 triangles)
 *
 * Creates a box mesh that bounds the input mesh.
 * Useful for very fast, coarse collision detection.
 *
 * @param mesh      Input mesh
 * @param use_obb   Use OBB (1) or AABB (0)
 * @param box_out   Output box mesh (must be initialized)
 * @return          GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_generate_box_proxy(const geomlib_trimesh *mesh, int use_obb, geomlib_trimesh *box_out);

/* ============================================================================
 * Mesh Analysis for Collision
 * ============================================================================ */

/**
 * Compute mesh quality metrics
 */
typedef struct geomlib_mesh_quality {
    float min_edge_length;      /* Shortest edge */
    float max_edge_length;      /* Longest edge */
    float avg_edge_length;      /* Average edge length */
    float min_triangle_area;    /* Smallest triangle */
    float max_triangle_area;    /* Largest triangle */
    size_t degenerate_count;    /* Number of degenerate triangles */
} geomlib_mesh_quality;

/**
 * Analyze mesh quality
 *
 * @param mesh          Mesh to analyze
 * @param quality_out   Output quality metrics
 * @return              GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_analyze_mesh_quality(const geomlib_trimesh *mesh, geomlib_mesh_quality *quality_out);

/**
 * Estimate collision detection cost
 *
 * Returns estimated relative cost of collision detection
 * (higher = more expensive). Useful for comparing proxy meshes.
 *
 * @param mesh      Mesh to analyze
 * @return          Relative cost estimate
 */
float geomlib_estimate_collision_cost(const geomlib_trimesh *mesh);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_COLLISION_H */
