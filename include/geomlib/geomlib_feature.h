/*
 * geomlib - Feature Detection
 * Sharp edge, ridge, valley, and corner detection algorithms
 */

#ifndef GEOMLIB_FEATURE_H
#define GEOMLIB_FEATURE_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include "geomlib_topology.h"
#include "geomlib_curvature.h"
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Feature Types and Structures
 * ============================================================================ */

/* Feature type classification */
typedef enum geomlib_feature_type {
    GEOMLIB_FEATURE_NONE = 0,
    GEOMLIB_FEATURE_SHARP_EDGE,     /* Dihedral angle exceeds threshold */
    GEOMLIB_FEATURE_RIDGE,          /* Positive curvature extremum */
    GEOMLIB_FEATURE_VALLEY,         /* Negative curvature extremum */
    GEOMLIB_FEATURE_CORNER,         /* Multiple sharp edges meet */
    GEOMLIB_FEATURE_BOUNDARY        /* Mesh boundary edge */
} geomlib_feature_type;

/* Feature edge descriptor */
typedef struct geomlib_feature_edge {
    uint32_t v0, v1;                /* Edge vertices */
    geomlib_feature_type type;
    float sharpness;                /* Dihedral angle (degrees) or curvature */
    uint32_t face0, face1;          /* Adjacent faces (face1 = UINT32_MAX if boundary) */
} geomlib_feature_edge;

/* Feature vertex descriptor (corner) */
typedef struct geomlib_feature_vertex {
    uint32_t vertex_index;
    geomlib_feature_type type;
    float sharpness;                /* Maximum dihedral angle at vertex */
    uint32_t edge_count;            /* Number of feature edges meeting here */
} geomlib_feature_vertex;

/* Complete feature set */
typedef struct geomlib_features {
    geomlib_feature_edge *edges;
    size_t edge_count;

    geomlib_feature_vertex *vertices;
    size_t vertex_count;

    /* Feature curves (chains of connected feature edges) */
    uint32_t **curves;              /* Array of vertex index arrays */
    size_t *curve_lengths;
    size_t curve_count;
} geomlib_features;

/* Feature detection parameters */
typedef struct geomlib_feature_params {
    /* Sharp edge detection */
    float sharp_angle_threshold;    /* Dihedral angle threshold (degrees, typical: 30-60) */

    /* Ridge/valley detection */
    bool detect_ridges_valleys;
    float curvature_threshold;       /* Minimum curvature for ridge/valley */

    /* Corner detection */
    bool detect_corners;
    uint32_t min_corner_edges;       /* Minimum sharp edges for corner (typical: 3) */

    /* Boundary handling */
    bool mark_boundaries;            /* Mark boundary edges as features */

    /* Filtering */
    float min_edge_length;           /* Ignore edges shorter than this */
    bool filter_small_features;      /* Remove isolated small features */
} geomlib_feature_params;

/* ============================================================================
 * Main Feature Detection
 * ============================================================================ */

/* Default detection parameters */
geomlib_feature_params geomlib_feature_params_default(void);

/* Detect all features in mesh */
geomlib_features geomlib_detect_features(
    const geomlib_trimesh *mesh,
    const geomlib_feature_params *params
);

/* Free features structure */
void geomlib_features_free(geomlib_features *features);

/* ============================================================================
 * Sharp Edge Detection
 * ============================================================================ */

/* Detect sharp edges by dihedral angle */
geomlib_feature_edge* geomlib_detect_sharp_edges(
    const geomlib_trimesh *mesh,
    float angle_threshold_deg,
    size_t *edge_count
);

/* Compute dihedral angle for edge */
float geomlib_compute_dihedral_angle(
    const geomlib_trimesh *mesh,
    uint32_t v0,
    uint32_t v1,
    uint32_t face0,
    uint32_t face1
);

/* Check if edge is sharp */
bool geomlib_is_edge_sharp(
    const geomlib_trimesh *mesh,
    uint32_t v0,
    uint32_t v1,
    float angle_threshold_deg
);

/* ============================================================================
 * Ridge and Valley Detection
 * ============================================================================ */

/* Detect ridges and valleys using curvature */
geomlib_feature_edge* geomlib_detect_ridges_valleys(
    const geomlib_trimesh *mesh,
    const geomlib_mesh_curvature *curvature,
    float threshold,
    size_t *edge_count
);

/* Check if edge is a ridge (high positive curvature) */
bool geomlib_is_ridge(
    const geomlib_trimesh *mesh,
    const geomlib_mesh_curvature *curvature,
    uint32_t v0,
    uint32_t v1,
    float threshold
);

/* Check if edge is a valley (high negative curvature) */
bool geomlib_is_valley(
    const geomlib_trimesh *mesh,
    const geomlib_mesh_curvature *curvature,
    uint32_t v0,
    uint32_t v1,
    float threshold
);

/* ============================================================================
 * Corner Detection
 * ============================================================================ */

/* Detect corner vertices (3+ feature edges meet) */
geomlib_feature_vertex* geomlib_detect_corners(
    const geomlib_trimesh *mesh,
    const geomlib_feature_edge *feature_edges,
    size_t edge_count,
    uint32_t min_edges,
    size_t *corner_count
);

/* Check if vertex is a corner */
bool geomlib_is_corner(
    const geomlib_feature_edge *feature_edges,
    size_t edge_count,
    uint32_t vertex_index,
    uint32_t min_edges
);

/* ============================================================================
 * Feature Curve Extraction
 * ============================================================================ */

/* Extract feature curves (chains of connected feature edges) */
int geomlib_extract_feature_curves(
    const geomlib_feature_edge *edges,
    size_t edge_count,
    uint32_t ***curves,              /* Output: array of vertex arrays */
    size_t **curve_lengths,          /* Output: length of each curve */
    size_t *curve_count              /* Output: number of curves */
);

/* Trace single feature curve starting from vertex */
int geomlib_trace_feature_curve(
    const geomlib_feature_edge *edges,
    size_t edge_count,
    uint32_t start_vertex,
    bool *visited_edges,
    uint32_t *curve,                 /* Output: vertex indices */
    size_t *curve_length
);

/* ============================================================================
 * Boundary Feature Detection
 * ============================================================================ */

/* Mark boundary edges as features */
geomlib_feature_edge* geomlib_mark_boundary_features(
    const geomlib_trimesh *mesh,
    size_t *edge_count
);

/* ============================================================================
 * Planar Region Segmentation
 * ============================================================================ */

/* Region descriptor */
typedef struct geomlib_planar_region {
    uint32_t *triangles;            /* Triangle indices in region */
    size_t triangle_count;

    geomlib_vec3 plane_origin;
    geomlib_vec3 plane_normal;

    float planarity;                /* How planar is this region [0,1] */
    float area;
} geomlib_planar_region;

/* Planar region segmentation result */
typedef struct geomlib_regions {
    geomlib_planar_region *regions;
    size_t region_count;

    uint32_t *triangle_labels;      /* Region ID for each triangle */
} geomlib_regions;

/* Segment mesh into planar regions */
geomlib_regions geomlib_segment_planar_regions(
    const geomlib_trimesh *mesh,
    float angle_threshold_deg,
    float distance_threshold
);

/* Free regions structure */
void geomlib_regions_free(geomlib_regions *regions);

/* ============================================================================
 * Feature Utilities
 * ============================================================================ */

/* Count feature edges at vertex */
uint32_t geomlib_count_feature_edges_at_vertex(
    const geomlib_feature_edge *edges,
    size_t edge_count,
    uint32_t vertex_index
);

/* Find all feature edges at vertex */
void geomlib_get_feature_edges_at_vertex(
    const geomlib_feature_edge *edges,
    size_t edge_count,
    uint32_t vertex_index,
    geomlib_feature_edge **vertex_edges,    /* Output */
    size_t *vertex_edge_count               /* Output */
);

/* Compute feature statistics */
typedef struct geomlib_feature_stats {
    size_t sharp_edge_count;
    size_t ridge_count;
    size_t valley_count;
    size_t corner_count;
    size_t boundary_edge_count;

    float avg_sharpness;
    float max_sharpness;
    float min_sharpness;
} geomlib_feature_stats;

geomlib_feature_stats geomlib_compute_feature_stats(
    const geomlib_features *features
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_FEATURE_H */
