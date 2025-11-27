/*
 * geomlib - Mesh Quality Analysis Implementation
 */

#include "geomlib/geomlib_quality.h"
#include "geomlib/geomlib_healing.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Default quality thresholds */
geomlib_quality_thresholds geomlib_quality_thresholds_default(void) {
    geomlib_quality_thresholds t;
    t.min_angle_threshold = 5.0f;      /* 5 degrees */
    t.max_angle_threshold = 170.0f;    /* 170 degrees */
    t.aspect_ratio_threshold = 10.0f;
    t.min_area_threshold = 1e-10f;
    t.sliver_threshold = 100.0f;
    return t;
}

/* Compute quality for a single triangle */
geomlib_triangle_quality geomlib_compute_triangle_quality(
    geomlib_vec3 p0,
    geomlib_vec3 p1,
    geomlib_vec3 p2,
    const geomlib_quality_thresholds *thresholds
) {
    geomlib_triangle_quality quality;
    memset(&quality, 0, sizeof(quality));

    /* Compute edges */
    geomlib_vec3 e0 = geomlib_vec3_sub(p1, p0);
    geomlib_vec3 e1 = geomlib_vec3_sub(p2, p1);
    geomlib_vec3 e2 = geomlib_vec3_sub(p0, p2);

    float len0 = geomlib_vec3_length(e0);
    float len1 = geomlib_vec3_length(e1);
    float len2 = geomlib_vec3_length(e2);

    /* Compute area */
    geomlib_vec3 cross = geomlib_vec3_cross(e0, geomlib_vec3_negate(e2));
    float area = geomlib_vec3_length(cross) * 0.5f;
    quality.area = area;

    /* Check for degenerate triangle */
    if (area < thresholds->min_area_threshold) {
        quality.is_degenerate = true;
        quality.shape_measure = 0.0f;
        quality.aspect_ratio = FLT_MAX;
        quality.radius_ratio = FLT_MAX;
        quality.min_angle = 0.0f;
        quality.max_angle = 180.0f;
        return quality;
    }

    /* Aspect ratio (longest edge / shortest edge) */
    float max_len = fmaxf(fmaxf(len0, len1), len2);
    float min_len = fminf(fminf(len0, len1), len2);
    quality.aspect_ratio = (min_len > 1e-10f) ? (max_len / min_len) : FLT_MAX;

    /* Check for sliver */
    if (quality.aspect_ratio > thresholds->sliver_threshold) {
        quality.is_sliver = true;
    }

    /* Compute angles using law of cosines */
    float angle0, angle1, angle2;

    if (len0 > 1e-10f && len2 > 1e-10f) {
        float cos_angle = -geomlib_vec3_dot(e0, e2) / (len0 * len2);
        cos_angle = fmaxf(-1.0f, fminf(1.0f, cos_angle));
        angle0 = acosf(cos_angle) * 180.0f / M_PI;
    } else {
        angle0 = 0.0f;
    }

    if (len0 > 1e-10f && len1 > 1e-10f) {
        float cos_angle = -geomlib_vec3_dot(e0, e1) / (len0 * len1);
        cos_angle = fmaxf(-1.0f, fminf(1.0f, cos_angle));
        angle1 = acosf(cos_angle) * 180.0f / M_PI;
    } else {
        angle1 = 0.0f;
    }

    angle2 = 180.0f - angle0 - angle1;

    quality.min_angle = fminf(fminf(angle0, angle1), angle2);
    quality.max_angle = fmaxf(fmaxf(angle0, angle1), angle2);

    /* Circumradius and inradius */
    float perimeter = len0 + len1 + len2;
    float inradius = (2.0f * area) / perimeter;

    /* Circumradius = (a*b*c) / (4*area) */
    float circumradius = (len0 * len1 * len2) / (4.0f * area);

    quality.radius_ratio = (inradius > 1e-10f) ? (circumradius / inradius) : FLT_MAX;

    /* Shape measure (1 = equilateral, 0 = degenerate) */
    /* Based on ratio of inradius to circumradius */
    /* For equilateral triangle: inradius/circumradius = 0.5 */
    if (circumradius > 1e-10f) {
        float ratio = inradius / circumradius;
        quality.shape_measure = ratio * 2.0f;  /* Normalize so equilateral = 1 */
    } else {
        quality.shape_measure = 0.0f;
    }

    return quality;
}

/* Compute quality statistics for entire mesh */
int geomlib_compute_mesh_quality(
    const geomlib_trimesh *mesh,
    const geomlib_quality_thresholds *thresholds,
    geomlib_mesh_quality_stats *stats
) {
    if (!mesh || !stats) return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_quality_thresholds default_thresholds = geomlib_quality_thresholds_default();
    if (!thresholds) thresholds = &default_thresholds;

    memset(stats, 0, sizeof(geomlib_mesh_quality_stats));

    /* Initialize min/max values */
    stats->min_triangle_area = FLT_MAX;
    stats->max_triangle_area = -FLT_MAX;
    stats->min_angle_global = FLT_MAX;
    stats->max_angle_global = -FLT_MAX;
    stats->min_aspect_ratio = FLT_MAX;
    stats->max_aspect_ratio = -FLT_MAX;
    stats->min_shape_measure = FLT_MAX;
    stats->min_edge_length = FLT_MAX;
    stats->max_edge_length = -FLT_MAX;
    stats->min_valence = INT32_MAX;
    stats->max_valence = 0;

    /* Compute per-triangle quality */
    float sum_area = 0.0f;
    float sum_min_angle = 0.0f;
    float sum_aspect_ratio = 0.0f;
    float sum_shape_measure = 0.0f;
    float sum_edge_length = 0.0f;
    float sum_edge_length_sq = 0.0f;
    int edge_count = 0;

    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        geomlib_vec3 p0 = mesh->vertices[i0];
        geomlib_vec3 p1 = mesh->vertices[i1];
        geomlib_vec3 p2 = mesh->vertices[i2];

        geomlib_triangle_quality q = geomlib_compute_triangle_quality(p0, p1, p2, thresholds);

        /* Update statistics */
        if (q.area < stats->min_triangle_area) stats->min_triangle_area = q.area;
        if (q.area > stats->max_triangle_area) stats->max_triangle_area = q.area;
        sum_area += q.area;

        if (q.min_angle < stats->min_angle_global) stats->min_angle_global = q.min_angle;
        if (q.max_angle > stats->max_angle_global) stats->max_angle_global = q.max_angle;
        sum_min_angle += q.min_angle;

        if (q.aspect_ratio < stats->min_aspect_ratio) stats->min_aspect_ratio = q.aspect_ratio;
        if (q.aspect_ratio > stats->max_aspect_ratio) stats->max_aspect_ratio = q.aspect_ratio;
        sum_aspect_ratio += q.aspect_ratio;

        if (q.shape_measure < stats->min_shape_measure) stats->min_shape_measure = q.shape_measure;
        sum_shape_measure += q.shape_measure;

        if (q.is_degenerate) stats->degenerate_count++;
        if (q.is_sliver) stats->sliver_count++;

        /* Edge lengths */
        float len0 = geomlib_vec3_length(geomlib_vec3_sub(p1, p0));
        float len1 = geomlib_vec3_length(geomlib_vec3_sub(p2, p1));
        float len2 = geomlib_vec3_length(geomlib_vec3_sub(p0, p2));

        for (int i = 0; i < 3; i++) {
            float len = (i == 0) ? len0 : (i == 1) ? len1 : len2;
            if (len < stats->min_edge_length) stats->min_edge_length = len;
            if (len > stats->max_edge_length) stats->max_edge_length = len;
            sum_edge_length += len;
            sum_edge_length_sq += len * len;
            edge_count++;
        }
    }

    /* Compute means */
    if (mesh->triangle_count > 0) {
        stats->mean_triangle_area = sum_area / mesh->triangle_count;
        stats->mean_min_angle = sum_min_angle / mesh->triangle_count;
        stats->mean_aspect_ratio = sum_aspect_ratio / mesh->triangle_count;
        stats->mean_shape_measure = sum_shape_measure / mesh->triangle_count;
    }

    if (edge_count > 0) {
        stats->mean_edge_length = sum_edge_length / edge_count;
        float variance = (sum_edge_length_sq / edge_count) -
                        (stats->mean_edge_length * stats->mean_edge_length);
        stats->edge_length_std_dev = sqrtf(fmaxf(0.0f, variance));
    }

    /* Build half-edge structure for topology analysis */
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, mesh);
    if (result == GEOMLIB_MESH_SUCCESS) {
        /* Count boundary and non-manifold edges */
        for (size_t he = 0; he < hemesh.halfedge_count; he++) {
            if (hemesh.halfedges[he].twin == UINT32_MAX) {
                stats->boundary_edges++;
            }
        }

        /* Compute vertex valence statistics */
        int *valences = calloc(mesh->vertex_count, sizeof(int));

        for (size_t he = 0; he < hemesh.halfedge_count; he++) {
            uint32_t v = hemesh.halfedges[he].vertex;
            if (v < mesh->vertex_count) {
                valences[v]++;
            }
        }

        float sum_valence = 0.0f;
        for (size_t v = 0; v < mesh->vertex_count; v++) {
            int val = valences[v];
            if (val < stats->min_valence) stats->min_valence = val;
            if (val > stats->max_valence) stats->max_valence = val;
            sum_valence += val;

            if (val != 6) {  /* Irregular for triangular mesh */
                stats->irregular_vertices++;
            }
        }

        stats->mean_valence = sum_valence / mesh->vertex_count;

        free(valences);
        geomlib_halfedge_free(&hemesh);
    }

    /* Compute overall quality score [0,1] */
    /* Based on shape measure and lack of bad triangles */
    float shape_score = stats->mean_shape_measure;
    float angle_score = (stats->min_angle_global > 10.0f) ? 1.0f : stats->min_angle_global / 10.0f;
    float degenerate_penalty = 1.0f - ((float)stats->degenerate_count / mesh->triangle_count);
    float sliver_penalty = 1.0f - ((float)stats->sliver_count / mesh->triangle_count);

    stats->quality_score = (shape_score + angle_score + degenerate_penalty + sliver_penalty) / 4.0f;

    return GEOMLIB_MESH_SUCCESS;
}

/* Get per-triangle quality array */
int geomlib_compute_per_triangle_quality(
    const geomlib_trimesh *mesh,
    const geomlib_quality_thresholds *thresholds,
    geomlib_triangle_quality **quality_array
) {
    if (!mesh || !quality_array) return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_quality_thresholds default_thresholds = geomlib_quality_thresholds_default();
    if (!thresholds) thresholds = &default_thresholds;

    *quality_array = malloc(mesh->triangle_count * sizeof(geomlib_triangle_quality));
    if (!*quality_array) return GEOMLIB_MESH_ERROR_MEMORY;

    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        (*quality_array)[t] = geomlib_compute_triangle_quality(
            mesh->vertices[i0],
            mesh->vertices[i1],
            mesh->vertices[i2],
            thresholds
        );
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* Identify problematic triangles */
int geomlib_identify_bad_triangles(
    const geomlib_trimesh *mesh,
    const geomlib_quality_thresholds *thresholds,
    uint32_t **bad_triangle_indices,
    size_t *bad_count
) {
    if (!mesh || !bad_triangle_indices || !bad_count) return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_quality_thresholds default_thresholds = geomlib_quality_thresholds_default();
    if (!thresholds) thresholds = &default_thresholds;

    /* First pass: count bad triangles */
    size_t count = 0;
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        geomlib_triangle_quality q = geomlib_compute_triangle_quality(
            mesh->vertices[i0], mesh->vertices[i1], mesh->vertices[i2], thresholds);

        if (q.is_degenerate || q.is_sliver ||
            q.min_angle < thresholds->min_angle_threshold ||
            q.max_angle > thresholds->max_angle_threshold) {
            count++;
        }
    }

    if (count == 0) {
        *bad_triangle_indices = NULL;
        *bad_count = 0;
        return GEOMLIB_MESH_SUCCESS;
    }

    /* Second pass: collect bad triangle indices */
    *bad_triangle_indices = malloc(count * sizeof(uint32_t));
    if (!*bad_triangle_indices) return GEOMLIB_MESH_ERROR_MEMORY;

    size_t idx = 0;
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        geomlib_triangle_quality q = geomlib_compute_triangle_quality(
            mesh->vertices[i0], mesh->vertices[i1], mesh->vertices[i2], thresholds);

        if (q.is_degenerate || q.is_sliver ||
            q.min_angle < thresholds->min_angle_threshold ||
            q.max_angle > thresholds->max_angle_threshold) {
            (*bad_triangle_indices)[idx++] = t;
        }
    }

    *bad_count = count;
    return GEOMLIB_MESH_SUCCESS;
}

/* Check if mesh meets quality criteria */
bool geomlib_check_mesh_quality(
    const geomlib_trimesh *mesh,
    const geomlib_quality_thresholds *thresholds
) {
    if (!mesh) return false;

    geomlib_mesh_quality_stats stats;
    if (geomlib_compute_mesh_quality(mesh, thresholds, &stats) != GEOMLIB_MESH_SUCCESS) {
        return false;
    }

    /* Pass if no degenerate triangles and quality score is reasonable */
    return (stats.degenerate_count == 0 && stats.quality_score > 0.5f);
}

/* Compute quality histogram */
int geomlib_compute_quality_histogram(
    const geomlib_trimesh *mesh,
    int bin_count,
    geomlib_quality_histogram *histogram
) {
    if (!mesh || !histogram || bin_count <= 0) return GEOMLIB_MESH_ERROR_INVALID;

    histogram->bin_count = bin_count;
    histogram->angle_bins = calloc(bin_count, sizeof(float));
    histogram->aspect_ratio_bins = calloc(bin_count, sizeof(float));

    if (!histogram->angle_bins || !histogram->aspect_ratio_bins) {
        geomlib_quality_histogram_free(histogram);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    geomlib_quality_thresholds thresholds = geomlib_quality_thresholds_default();

    /* Collect data */
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        geomlib_triangle_quality q = geomlib_compute_triangle_quality(
            mesh->vertices[i0], mesh->vertices[i1], mesh->vertices[i2], &thresholds);

        /* Angle histogram (0-180 degrees) */
        int angle_bin = (int)((q.min_angle / 180.0f) * bin_count);
        if (angle_bin >= 0 && angle_bin < bin_count) {
            histogram->angle_bins[angle_bin]++;
        }

        /* Aspect ratio histogram (log scale, 1-100) */
        float log_ar = log10f(fminf(q.aspect_ratio, 100.0f));
        int ar_bin = (int)((log_ar / 2.0f) * bin_count);
        if (ar_bin >= 0 && ar_bin < bin_count) {
            histogram->aspect_ratio_bins[ar_bin]++;
        }
    }

    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_quality_histogram_free(geomlib_quality_histogram *histogram) {
    if (!histogram) return;
    if (histogram->angle_bins) {
        free(histogram->angle_bins);
        histogram->angle_bins = NULL;
    }
    if (histogram->aspect_ratio_bins) {
        free(histogram->aspect_ratio_bins);
        histogram->aspect_ratio_bins = NULL;
    }
    histogram->bin_count = 0;
}
