/*
 * geomlib - Mesh Quality Analysis
 * Metrics for evaluating triangle mesh quality
 */

#ifndef GEOMLIB_QUALITY_H
#define GEOMLIB_QUALITY_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Triangle quality metrics */
typedef struct geomlib_triangle_quality {
    float aspect_ratio;        /* Edge length ratio (max/min) */
    float min_angle;           /* Minimum interior angle (degrees) */
    float max_angle;           /* Maximum interior angle (degrees) */
    float area;                /* Triangle area */
    float radius_ratio;        /* Circumradius / inradius */
    float shape_measure;       /* Quality metric in [0,1], 1=equilateral */
    bool is_degenerate;        /* True if area is near zero */
    bool is_sliver;            /* True if very flat (aspect ratio > threshold) */
} geomlib_triangle_quality;

/* Mesh quality statistics */
typedef struct geomlib_mesh_quality_stats {
    /* Triangle statistics */
    float min_triangle_area;
    float max_triangle_area;
    float mean_triangle_area;

    float min_angle_global;     /* Worst angle in entire mesh */
    float max_angle_global;
    float mean_min_angle;       /* Average of minimum angles */

    float min_aspect_ratio;
    float max_aspect_ratio;
    float mean_aspect_ratio;

    float min_shape_measure;
    float mean_shape_measure;

    /* Edge statistics */
    float min_edge_length;
    float max_edge_length;
    float mean_edge_length;
    float edge_length_std_dev;

    /* Mesh issues */
    uint32_t degenerate_count;  /* Triangles with near-zero area */
    uint32_t sliver_count;      /* Very flat triangles */
    uint32_t non_manifold_edges; /* Edges shared by >2 triangles */
    uint32_t boundary_edges;    /* Edges with only 1 triangle */

    /* Vertex valence statistics */
    int min_valence;            /* Minimum vertex degree */
    int max_valence;            /* Maximum vertex degree */
    float mean_valence;
    uint32_t irregular_vertices; /* Vertices with valence != 6 (for triangular meshes) */

    /* Overall quality score [0,1] */
    float quality_score;        /* Composite quality metric */
} geomlib_mesh_quality_stats;

/* Quality thresholds */
typedef struct geomlib_quality_thresholds {
    float min_angle_threshold;      /* Degrees, default: 5 */
    float max_angle_threshold;      /* Degrees, default: 170 */
    float aspect_ratio_threshold;   /* Default: 10 */
    float min_area_threshold;       /* Default: 1e-10 */
    float sliver_threshold;         /* Aspect ratio, default: 100 */
} geomlib_quality_thresholds;

/* Default quality thresholds */
geomlib_quality_thresholds geomlib_quality_thresholds_default(void);

/* Compute quality for a single triangle */
geomlib_triangle_quality geomlib_compute_triangle_quality(
    geomlib_vec3 p0,
    geomlib_vec3 p1,
    geomlib_vec3 p2,
    const geomlib_quality_thresholds *thresholds
);

/* Compute quality statistics for entire mesh */
int geomlib_compute_mesh_quality(
    const geomlib_trimesh *mesh,
    const geomlib_quality_thresholds *thresholds,
    geomlib_mesh_quality_stats *stats
);

/* Get per-triangle quality array */
int geomlib_compute_per_triangle_quality(
    const geomlib_trimesh *mesh,
    const geomlib_quality_thresholds *thresholds,
    geomlib_triangle_quality **quality_array  /* Output: allocated array */
);

/* Identify problematic triangles */
int geomlib_identify_bad_triangles(
    const geomlib_trimesh *mesh,
    const geomlib_quality_thresholds *thresholds,
    uint32_t **bad_triangle_indices,  /* Output: indices of bad triangles */
    size_t *bad_count                 /* Output: number of bad triangles */
);

/* Check if mesh meets quality criteria */
bool geomlib_check_mesh_quality(
    const geomlib_trimesh *mesh,
    const geomlib_quality_thresholds *thresholds
);

/* Compute mesh histogram data for visualization */
typedef struct geomlib_quality_histogram {
    float *angle_bins;        /* Histogram of angles */
    float *aspect_ratio_bins; /* Histogram of aspect ratios */
    int bin_count;            /* Number of bins */
} geomlib_quality_histogram;

int geomlib_compute_quality_histogram(
    const geomlib_trimesh *mesh,
    int bin_count,
    geomlib_quality_histogram *histogram
);

void geomlib_quality_histogram_free(geomlib_quality_histogram *histogram);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_QUALITY_H */
