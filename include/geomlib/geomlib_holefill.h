/*
 * geomlib - Hole Filling and Mesh Repair
 * Automatic hole detection and filling algorithms
 */

#ifndef GEOMLIB_HOLEFILL_H
#define GEOMLIB_HOLEFILL_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include "geomlib_topology.h"
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Hole Detection and Classification
 * ============================================================================ */

/* Hole descriptor */
typedef struct geomlib_hole {
    uint32_t *boundary_vertices;    /* Ordered vertices around hole */
    size_t vertex_count;

    geomlib_vec3 centroid;          /* Center of hole */
    geomlib_vec3 normal;            /* Average normal at boundary */
    float perimeter;                /* Total boundary length */
    float area_estimate;            /* Estimated hole area */

    bool is_planar;                 /* Is hole roughly planar? */
    float planarity;                /* Planarity measure [0,1] */
} geomlib_hole;

/* Hole detection result */
typedef struct geomlib_holes {
    geomlib_hole *holes;
    size_t hole_count;
} geomlib_holes;

/* Detect all holes in mesh */
geomlib_holes geomlib_detect_holes(const geomlib_trimesh *mesh);

/* Free holes structure */
void geomlib_holes_free(geomlib_holes *holes);

/* Classify single hole */
geomlib_hole geomlib_classify_hole(
    const geomlib_trimesh *mesh,
    const geomlib_boundary_loop *loop
);

/* ============================================================================
 * Hole Filling Methods
 * ============================================================================ */

/* Hole filling method */
typedef enum geomlib_fill_method {
    GEOMLIB_FILL_ADVANCING_FRONT,   /* Advancing front triangulation */
    GEOMLIB_FILL_MINIMAL_SURFACE,   /* Minimal surface with smoothing */
    GEOMLIB_FILL_TANGENT_PLANE,     /* Flat fill in tangent plane */
    GEOMLIB_FILL_NATURAL             /* Feature-preserving natural fill */
} geomlib_fill_method;

/* Hole filling parameters */
typedef struct geomlib_fill_params {
    geomlib_fill_method method;

    /* Smoothing parameters */
    uint32_t smooth_iterations;     /* Number of smoothing passes */
    float smooth_lambda;            /* Smoothing factor [0,1] */

    /* Quality parameters */
    float max_edge_length;          /* Maximum edge length (0 = auto) */
    float min_angle_deg;            /* Minimum triangle angle */

    /* Feature preservation */
    bool preserve_boundary_edges;   /* Don't subdivide boundary */
    bool fair_interior;             /* Apply fairing to interior */
} geomlib_fill_params;

/* Default fill parameters */
geomlib_fill_params geomlib_fill_params_default(geomlib_fill_method method);

/* Fill single hole */
int geomlib_fill_hole(
    geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    const geomlib_fill_params *params
);

/* Fill all holes in mesh */
int geomlib_fill_all_holes(
    geomlib_trimesh *mesh,
    const geomlib_fill_params *params
);

/* ============================================================================
 * Advancing Front Triangulation
 * ============================================================================ */

/* Advancing front context */
typedef struct geomlib_advancing_front {
    uint32_t *front_edges;          /* Current front edges (v0, v1 pairs) */
    size_t front_count;

    uint32_t *interior_vertices;    /* Newly created vertices */
    size_t interior_vertex_count;

    uint32_t *new_triangles;        /* Generated triangle indices */
    size_t triangle_count;
} geomlib_advancing_front;

/* Fill hole using advancing front */
int geomlib_fill_advancing_front(
    geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    const geomlib_fill_params *params
);

/* ============================================================================
 * Minimal Surface Filling
 * ============================================================================ */

/* Fill hole with minimal surface (smoothed) */
int geomlib_fill_minimal_surface(
    geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    const geomlib_fill_params *params
);

/* ============================================================================
 * Tangent Plane Filling
 * ============================================================================ */

/* Project hole to plane and triangulate */
int geomlib_fill_tangent_plane(
    geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    const geomlib_fill_params *params
);

/* ============================================================================
 * Smoothing and Fairing
 * ============================================================================ */

/* Fair (smooth) filled region */
int geomlib_fair_filled_region(
    geomlib_trimesh *mesh,
    const uint32_t *region_vertices,
    size_t vertex_count,
    uint32_t iterations,
    float lambda
);

/* Laplacian smoothing for hole interior */
int geomlib_smooth_hole_interior(
    geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    const uint32_t *interior_vertices,
    size_t interior_count,
    uint32_t iterations
);

/* ============================================================================
 * Hole Utilities
 * ============================================================================ */

/* Compute hole centroid */
geomlib_vec3 geomlib_hole_centroid(
    const geomlib_trimesh *mesh,
    const geomlib_hole *hole
);

/* Compute hole normal (average of boundary normals) */
geomlib_vec3 geomlib_hole_normal(
    const geomlib_trimesh *mesh,
    const geomlib_hole *hole
);

/* Estimate hole area */
float geomlib_hole_area_estimate(
    const geomlib_trimesh *mesh,
    const geomlib_hole *hole
);

/* Check if hole is approximately planar */
bool geomlib_is_hole_planar(
    const geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    float threshold
);

/* Compute best-fit plane for hole */
void geomlib_hole_best_fit_plane(
    const geomlib_trimesh *mesh,
    const geomlib_hole *hole,
    geomlib_vec3 *plane_origin,
    geomlib_vec3 *plane_normal
);

/* ============================================================================
 * Polygon Triangulation (2D)
 * ============================================================================ */

/* Triangulate 2D polygon (for planar holes) */
int geomlib_triangulate_polygon_2d(
    const geomlib_vec2 *points,
    size_t point_count,
    uint32_t *triangles,  /* Output: 3*N indices */
    size_t *triangle_count
);

/* Ear clipping triangulation */
int geomlib_ear_clipping(
    const geomlib_vec2 *points,
    size_t point_count,
    uint32_t *triangles,
    size_t *triangle_count
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_HOLEFILL_H */
