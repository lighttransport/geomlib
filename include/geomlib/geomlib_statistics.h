/*
 * geomlib - Mesh Statistics and Validation
 * Comprehensive mesh analysis and diagnostic tools
 */

#ifndef GEOMLIB_STATISTICS_H
#define GEOMLIB_STATISTICS_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Edge Statistics
 * ============================================================================ */

typedef struct geomlib_edge_stats {
    float min_length;
    float max_length;
    float avg_length;
    float median_length;
    float std_dev;

    /* Histogram */
    uint32_t *histogram;
    size_t bin_count;
    float bin_width;
} geomlib_edge_stats;

geomlib_edge_stats geomlib_compute_edge_stats(
    const geomlib_trimesh *mesh,
    size_t histogram_bins
);

void geomlib_edge_stats_free(geomlib_edge_stats *stats);

/* ============================================================================
 * Angle Statistics
 * ============================================================================ */

typedef struct geomlib_angle_stats {
    float min_angle_deg;
    float max_angle_deg;
    float avg_angle_deg;

    /* Histogram */
    uint32_t *histogram;
    size_t bin_count;
    float bin_width;

    /* Quality metrics */
    uint32_t degenerate_count;      /* Angles < 1 degree */
    uint32_t poor_quality_count;    /* Angles < 10 or > 170 degrees */
} geomlib_angle_stats;

geomlib_angle_stats geomlib_compute_angle_stats(
    const geomlib_trimesh *mesh,
    size_t histogram_bins
);

void geomlib_angle_stats_free(geomlib_angle_stats *stats);

/* ============================================================================
 * Triangle Area Statistics
 * ============================================================================ */

typedef struct geomlib_area_stats {
    float min_area;
    float max_area;
    float avg_area;
    float total_area;
    float std_dev;

    /* Histogram */
    uint32_t *histogram;
    size_t bin_count;
    float bin_width;
} geomlib_area_stats;

geomlib_area_stats geomlib_compute_area_stats(
    const geomlib_trimesh *mesh,
    size_t histogram_bins
);

void geomlib_area_stats_free(geomlib_area_stats *stats);

/* ============================================================================
 * Aspect Ratio Statistics
 * ============================================================================ */

typedef struct geomlib_aspect_ratio_stats {
    float min_ratio;
    float max_ratio;
    float avg_ratio;

    /* Histogram */
    uint32_t *histogram;
    size_t bin_count;
    float bin_width;

    /* Quality counts */
    uint32_t good_count;            /* Ratio < 2 */
    uint32_t acceptable_count;      /* Ratio 2-5 */
    uint32_t poor_count;            /* Ratio > 5 */
} geomlib_aspect_ratio_stats;

geomlib_aspect_ratio_stats geomlib_compute_aspect_ratio_stats(
    const geomlib_trimesh *mesh,
    size_t histogram_bins
);

void geomlib_aspect_ratio_stats_free(geomlib_aspect_ratio_stats *stats);

/* ============================================================================
 * Comprehensive Mesh Statistics
 * ============================================================================ */

typedef struct geomlib_mesh_stats {
    /* Basic counts */
    size_t vertex_count;
    size_t triangle_count;
    size_t edge_count;

    /* Bounding box */
    geomlib_aabb bbox;
    geomlib_vec3 centroid;

    /* Geometric properties */
    float surface_area;
    float volume;               /* For closed meshes */
    bool is_closed;

    /* Edge statistics */
    geomlib_edge_stats edges;

    /* Angle statistics */
    geomlib_angle_stats angles;

    /* Area statistics */
    geomlib_area_stats areas;

    /* Aspect ratio statistics */
    geomlib_aspect_ratio_stats aspect_ratios;
} geomlib_mesh_stats;

geomlib_mesh_stats geomlib_compute_mesh_stats(
    const geomlib_trimesh *mesh,
    size_t histogram_bins
);

void geomlib_mesh_stats_free(geomlib_mesh_stats *stats);

/* ============================================================================
 * Mesh Validation
 * ============================================================================ */

typedef enum geomlib_validation_level {
    GEOMLIB_VALIDATE_BASIC,         /* Structure, indices, repeated indices, and zero area */
    GEOMLIB_VALIDATE_STANDARD,      /* Also check coordinates, unused vertices, and topology */
    GEOMLIB_VALIDATE_STRICT         /* Also check intersections and reject zero area */
} geomlib_validation_level;

typedef enum geomlib_validation_issue_type {
    GEOMLIB_ISSUE_DEGENERATE_TRIANGLE,
    GEOMLIB_ISSUE_DUPLICATE_VERTEX,
    GEOMLIB_ISSUE_UNUSED_VERTEX,
    GEOMLIB_ISSUE_INVALID_INDEX,
    GEOMLIB_ISSUE_NON_MANIFOLD_EDGE,
    GEOMLIB_ISSUE_NON_MANIFOLD_VERTEX,
    GEOMLIB_ISSUE_INCONSISTENT_ORIENTATION,
    GEOMLIB_ISSUE_SELF_INTERSECTION,
    GEOMLIB_ISSUE_ZERO_AREA_TRIANGLE,
    GEOMLIB_ISSUE_INVERTED_TRIANGLE,
    GEOMLIB_ISSUE_INVALID_MESH,
    GEOMLIB_ISSUE_NONFINITE_VERTEX
} geomlib_validation_issue_type;

typedef struct geomlib_validation_issue {
    geomlib_validation_issue_type type;
    uint32_t element_index;     /* Triangle, vertex, or edge index */
    char message[256];
} geomlib_validation_issue;

typedef struct geomlib_validation_result {
    bool is_valid;
    geomlib_validation_level level;

    geomlib_validation_issue *issues;
    size_t issue_count;

    /* Summary counts */
    uint32_t error_count;
    uint32_t warning_count;
} geomlib_validation_result;

geomlib_validation_result geomlib_validate_mesh(
    const geomlib_trimesh *mesh,
    geomlib_validation_level level
);

void geomlib_validation_result_free(geomlib_validation_result *result);

/* ============================================================================
 * Mesh Comparison
 * ============================================================================ */

typedef struct geomlib_mesh_difference {
    /* Vertex differences */
    size_t vertex_count_diff;
    size_t triangle_count_diff;

    /* Geometric differences */
    float hausdorff_distance;           /* Max bidirectional vertex-to-surface sample */
    float avg_distance;                 /* Average distance */
    float rms_distance;                 /* Root mean square */

    /* Volume difference */
    float volume_diff;
    float volume_diff_percent;
} geomlib_mesh_difference;

/* Distance statistics sample each mesh's vertices against the other mesh's triangles.
 * Two empty meshes have zero distance; a missing target surface yields FLT_MAX. */
geomlib_mesh_difference geomlib_compare_meshes(
    const geomlib_trimesh *mesh1,
    const geomlib_trimesh *mesh2
);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/* Compute triangle aspect ratio (longest edge / shortest edge) */
float geomlib_triangle_aspect_ratio(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2
);

/* Check if triangle is degenerate */
bool geomlib_is_degenerate_triangle(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    float threshold
);

/* Compute angle between three points (in degrees) */
float geomlib_compute_angle(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_STATISTICS_H */
