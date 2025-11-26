/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Mesh Healing and Repair
 * Algorithms for fixing mesh topology and geometry issues
 *
 * References:
 * - "Consistent Normal Orientation for Polygonal Meshes" (Borodin et al., 2004)
 * - "Robust Inside-Outside Segmentation using Generalized Winding Numbers" (Jacobson et al., 2013)
 * - "Surface Reconstruction from Unorganized Points" (Hoppe et al., SIGGRAPH 1992)
 * - "Mesh Repairing Using Topology Graphs" (Barequet & Kumar, 2005)
 */

#ifndef GEOMLIB_HEALING_H
#define GEOMLIB_HEALING_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include "geomlib_bvh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Half-Edge Data Structure
 * ============================================================================ */

/**
 * Half-edge structure for mesh topology
 *
 * Each edge is represented by two half-edges in opposite directions.
 * This allows efficient traversal of mesh connectivity.
 */
typedef struct geomlib_halfedge {
    uint32_t vertex;          /* Vertex this half-edge points to */
    uint32_t face;            /* Face this half-edge belongs to */
    uint32_t next;            /* Next half-edge in face (CCW) */
    uint32_t twin;            /* Opposite half-edge (or UINT32_MAX if boundary) */
} geomlib_halfedge;

/**
 * Half-edge mesh representation
 */
typedef struct geomlib_halfedge_mesh {
    geomlib_halfedge *halfedges;  /* Array of half-edges */
    size_t halfedge_count;        /* Number of half-edges */
    size_t halfedge_capacity;     /* Allocated capacity */

    uint32_t *vertex_halfedge;    /* First outgoing half-edge per vertex */
    uint32_t *face_halfedge;      /* First half-edge per face */

    const geomlib_trimesh *mesh;  /* Reference to original mesh */
} geomlib_halfedge_mesh;

/* ============================================================================
 * Half-Edge Operations
 * ============================================================================ */

/**
 * Build half-edge structure from triangle mesh
 *
 * @param hemesh    Half-edge mesh to initialize
 * @param mesh      Input triangle mesh
 * @return          GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_halfedge_build(
    geomlib_halfedge_mesh *hemesh,
    const geomlib_trimesh *mesh
);

/**
 * Free half-edge mesh resources
 *
 * @param hemesh    Half-edge mesh to free
 */
void geomlib_halfedge_free(geomlib_halfedge_mesh *hemesh);

/**
 * Get statistics about half-edge mesh
 *
 * @param hemesh                Half-edge mesh
 * @param boundary_edges_out    Output: number of boundary edges (can be NULL)
 * @param manifold_out          Output: true if mesh is manifold (can be NULL)
 * @param genus_out             Output: mesh genus (can be NULL)
 */
void geomlib_halfedge_get_stats(
    const geomlib_halfedge_mesh *hemesh,
    size_t *boundary_edges_out,
    bool *manifold_out,
    int *genus_out
);

/* ============================================================================
 * Face Orientation Fixing
 * ============================================================================ */

/**
 * Face orientation strategy
 */
typedef enum geomlib_orientation_strategy {
    GEOMLIB_ORIENT_PROPAGATION,     /* Propagate from consistent seed face */
    GEOMLIB_ORIENT_MAJORITY_VOTE,   /* Align to majority orientation */
    GEOMLIB_ORIENT_INSIDE_OUT,      /* Fix based on inside/outside detection */
    GEOMLIB_ORIENT_HYBRID           /* Combination of methods for robustness */
} geomlib_orientation_strategy;

/**
 * Face orientation parameters
 */
typedef struct geomlib_orientation_params {
    geomlib_orientation_strategy strategy;

    /* Inside/outside detection parameters */
    int ray_samples;              /* Number of ray samples (default 32) */
    float ray_offset;             /* Ray origin offset from surface (default 0.001) */

    /* Propagation parameters */
    bool flip_if_majority_inside; /* Flip all faces if majority point inward */
} geomlib_orientation_params;

/**
 * Create default orientation parameters
 *
 * @param strategy  Orientation strategy
 * @return          Default parameters for strategy
 */
geomlib_orientation_params geomlib_orientation_params_default(
    geomlib_orientation_strategy strategy
);

/**
 * Fix face orientations to be consistent
 *
 * Uses half-edge structure to propagate consistent winding order.
 * Optionally uses ray casting to determine inside/outside.
 *
 * @param mesh      Mesh to fix (will be modified in-place)
 * @param params    Orientation parameters
 * @param hemesh    Optional pre-built half-edge mesh (can be NULL)
 * @return          Number of faces flipped, or negative error code
 */
int geomlib_fix_face_orientation(
    geomlib_trimesh *mesh,
    const geomlib_orientation_params *params,
    geomlib_halfedge_mesh *hemesh
);

/**
 * Detect if face is pointing inward or outward
 *
 * Uses Monte Carlo ray casting to estimate if face normal
 * points toward interior or exterior of mesh.
 *
 * @param mesh          Mesh containing face
 * @param face_idx      Face index to test
 * @param bvh           Optional BVH for acceleration (can be NULL)
 * @param ray_samples   Number of ray samples
 * @return              Positive if outward, negative if inward, 0 if uncertain
 */
int geomlib_detect_face_orientation(
    const geomlib_trimesh *mesh,
    uint32_t face_idx,
    const geomlib_bvh *bvh,
    int ray_samples
);

/* ============================================================================
 * Mesh Healing Operations
 * ============================================================================ */

/**
 * Remove degenerate triangles
 *
 * Removes triangles with:
 * - Zero area
 * - Duplicate vertices
 * - Co-linear vertices
 *
 * @param mesh              Mesh to clean
 * @param area_threshold    Minimum triangle area (default 1e-10)
 * @return                  Number of triangles removed
 */
int geomlib_remove_degenerate_triangles(
    geomlib_trimesh *mesh,
    float area_threshold
);

/**
 * Remove duplicate vertices
 *
 * Merges vertices that are within distance threshold.
 * Updates triangle indices accordingly.
 *
 * @param mesh              Mesh to clean
 * @param distance_threshold Distance threshold for merging (default 1e-6)
 * @return                  Number of vertices merged
 */
int geomlib_remove_duplicate_vertices(
    geomlib_trimesh *mesh,
    float distance_threshold
);

/**
 * Fill holes in mesh
 *
 * Detects boundary loops and fills them with triangles.
 * Uses advancing front triangulation for hole filling.
 *
 * @param mesh              Mesh to repair
 * @param max_hole_edges    Maximum hole size to fill (default 100)
 * @return                  Number of holes filled
 */
int geomlib_fill_holes(
    geomlib_trimesh *mesh,
    int max_hole_edges
);

/**
 * Comprehensive mesh repair
 *
 * Applies multiple healing operations in optimal order:
 * 1. Remove duplicate vertices
 * 2. Remove degenerate triangles
 * 3. Fix face orientations
 * 4. Fill holes (optional)
 *
 * @param mesh              Mesh to repair
 * @param fill_holes        Whether to fill holes
 * @return                  GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_repair_mesh(
    geomlib_trimesh *mesh,
    bool fill_holes
);

/* ============================================================================
 * Mesh Quality Analysis
 * ============================================================================ */

/**
 * Mesh topology issues
 */
typedef struct geomlib_mesh_issues {
    size_t degenerate_triangles;   /* Number of degenerate triangles */
    size_t duplicate_vertices;     /* Number of duplicate vertices */
    size_t boundary_edges;         /* Number of boundary edges */
    size_t non_manifold_edges;     /* Number of non-manifold edges */
    size_t non_manifold_vertices;  /* Number of non-manifold vertices */
    size_t inconsistent_faces;     /* Faces with inconsistent orientation */
    bool is_manifold;              /* True if mesh is manifold */
    bool is_watertight;            /* True if mesh is closed (no boundaries) */
    bool is_oriented;              /* True if all faces consistently oriented */
} geomlib_mesh_issues;

/**
 * Analyze mesh for topology issues
 *
 * @param mesh      Mesh to analyze
 * @param issues    Output: detected issues
 * @return          GEOMLIB_MESH_SUCCESS on success
 */
int geomlib_analyze_mesh_topology(
    const geomlib_trimesh *mesh,
    geomlib_mesh_issues *issues
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_HEALING_H */
