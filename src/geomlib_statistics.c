/*
 * geomlib - Mesh Statistics Implementation (Basic)
 */

#include "geomlib/geomlib_statistics.h"
#include "geomlib/geomlib_topology.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

float geomlib_triangle_aspect_ratio(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2
) {
    float e0 = geomlib_vec3_length(geomlib_vec3_sub(v1, v0));
    float e1 = geomlib_vec3_length(geomlib_vec3_sub(v2, v1));
    float e2 = geomlib_vec3_length(geomlib_vec3_sub(v0, v2));

    float max_edge = fmaxf(e0, fmaxf(e1, e2));
    float min_edge = fminf(e0, fminf(e1, e2));

    return (min_edge > 1e-10f) ? (max_edge / min_edge) : FLT_MAX;
}

bool geomlib_is_degenerate_triangle(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    float threshold
) {
    geomlib_vec3 e0 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e1 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 n = geomlib_vec3_cross(e0, e1);
    float area = geomlib_vec3_length(n) * 0.5f;
    return area < threshold;
}

float geomlib_compute_angle(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2
) {
    geomlib_vec3 e1 = geomlib_vec3_normalize(geomlib_vec3_sub(v0, v1));
    geomlib_vec3 e2 = geomlib_vec3_normalize(geomlib_vec3_sub(v2, v1));
    float dot = geomlib_vec3_dot(e1, e2);
    dot = fmaxf(-1.0f, fminf(1.0f, dot));  /* Clamp */
    return acosf(dot) * 180.0f / 3.14159265359f;
}

/* ============================================================================
 * Edge Statistics
 * ============================================================================ */

geomlib_edge_stats geomlib_compute_edge_stats(
    const geomlib_trimesh *mesh,
    size_t histogram_bins
) {
    geomlib_edge_stats stats;
    stats.histogram = NULL;
    stats.bin_count = histogram_bins;

    /* Build edge mesh */
    geomlib_edge_mesh edge_mesh = geomlib_build_edge_mesh(mesh);

    float *lengths = (float*)malloc(edge_mesh.edge_count * sizeof(float));
    stats.min_length = FLT_MAX;
    stats.max_length = 0.0f;
    float sum = 0.0f;

    for (size_t i = 0; i < edge_mesh.edge_count; i++) {
        geomlib_vec3 v0 = mesh->vertices[edge_mesh.edges[i].v0];
        geomlib_vec3 v1 = mesh->vertices[edge_mesh.edges[i].v1];
        float len = geomlib_vec3_length(geomlib_vec3_sub(v1, v0));
        lengths[i] = len;
        sum += len;
        if (len < stats.min_length) stats.min_length = len;
        if (len > stats.max_length) stats.max_length = len;
    }

    stats.avg_length = sum / (float)edge_mesh.edge_count;

    /* Compute median and standard deviation */
    /* TODO: Sort for median */
    stats.median_length = stats.avg_length;

    float var_sum = 0.0f;
    for (size_t i = 0; i < edge_mesh.edge_count; i++) {
        float diff = lengths[i] - stats.avg_length;
        var_sum += diff * diff;
    }
    stats.std_dev = sqrtf(var_sum / (float)edge_mesh.edge_count);

    /* Build histogram */
    if (histogram_bins > 0) {
        stats.histogram = (uint32_t*)calloc(histogram_bins, sizeof(uint32_t));
        stats.bin_width = (stats.max_length - stats.min_length) / (float)histogram_bins;

        for (size_t i = 0; i < edge_mesh.edge_count; i++) {
            size_t bin = (size_t)((lengths[i] - stats.min_length) / stats.bin_width);
            if (bin >= histogram_bins) bin = histogram_bins - 1;
            stats.histogram[bin]++;
        }
    }

    free(lengths);
    geomlib_edge_mesh_free(&edge_mesh);

    return stats;
}

void geomlib_edge_stats_free(geomlib_edge_stats *stats) {
    free(stats->histogram);
}

/* ============================================================================
 * Angle Statistics
 * ============================================================================ */

geomlib_angle_stats geomlib_compute_angle_stats(
    const geomlib_trimesh *mesh,
    size_t histogram_bins
) {
    geomlib_angle_stats stats;
    stats.histogram = NULL;
    stats.bin_count = histogram_bins;
    stats.degenerate_count = 0;
    stats.poor_quality_count = 0;

    size_t angle_count = mesh->triangle_count * 3;
    float *angles = (float*)malloc(angle_count * sizeof(float));

    stats.min_angle_deg = 180.0f;
    stats.max_angle_deg = 0.0f;
    float sum = 0.0f;

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 v0 = mesh->vertices[mesh->indices[i * 3 + 0]];
        geomlib_vec3 v1 = mesh->vertices[mesh->indices[i * 3 + 1]];
        geomlib_vec3 v2 = mesh->vertices[mesh->indices[i * 3 + 2]];

        float a0 = geomlib_compute_angle(v1, v0, v2);
        float a1 = geomlib_compute_angle(v0, v1, v2);
        float a2 = geomlib_compute_angle(v0, v2, v1);

        angles[i * 3 + 0] = a0;
        angles[i * 3 + 1] = a1;
        angles[i * 3 + 2] = a2;

        for (int j = 0; j < 3; j++) {
            float a = angles[i * 3 + j];
            sum += a;
            if (a < stats.min_angle_deg) stats.min_angle_deg = a;
            if (a > stats.max_angle_deg) stats.max_angle_deg = a;

            if (a < 1.0f) stats.degenerate_count++;
            if (a < 10.0f || a > 170.0f) stats.poor_quality_count++;
        }
    }

    stats.avg_angle_deg = sum / (float)angle_count;

    /* Build histogram */
    if (histogram_bins > 0) {
        stats.histogram = (uint32_t*)calloc(histogram_bins, sizeof(uint32_t));
        stats.bin_width = 180.0f / (float)histogram_bins;

        for (size_t i = 0; i < angle_count; i++) {
            size_t bin = (size_t)(angles[i] / stats.bin_width);
            if (bin >= histogram_bins) bin = histogram_bins - 1;
            stats.histogram[bin]++;
        }
    }

    free(angles);
    return stats;
}

void geomlib_angle_stats_free(geomlib_angle_stats *stats) {
    free(stats->histogram);
}

/* ============================================================================
 * Area Statistics
 * ============================================================================ */

geomlib_area_stats geomlib_compute_area_stats(
    const geomlib_trimesh *mesh,
    size_t histogram_bins
) {
    geomlib_area_stats stats;
    stats.histogram = NULL;
    stats.bin_count = histogram_bins;

    float *areas = (float*)malloc(mesh->triangle_count * sizeof(float));
    stats.min_area = FLT_MAX;
    stats.max_area = 0.0f;
    stats.total_area = 0.0f;

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 v0 = mesh->vertices[mesh->indices[i * 3 + 0]];
        geomlib_vec3 v1 = mesh->vertices[mesh->indices[i * 3 + 1]];
        geomlib_vec3 v2 = mesh->vertices[mesh->indices[i * 3 + 2]];

        geomlib_vec3 e0 = geomlib_vec3_sub(v1, v0);
        geomlib_vec3 e1 = geomlib_vec3_sub(v2, v0);
        float area = geomlib_vec3_length(geomlib_vec3_cross(e0, e1)) * 0.5f;

        areas[i] = area;
        stats.total_area += area;
        if (area < stats.min_area) stats.min_area = area;
        if (area > stats.max_area) stats.max_area = area;
    }

    stats.avg_area = stats.total_area / (float)mesh->triangle_count;

    /* Standard deviation */
    float var_sum = 0.0f;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        float diff = areas[i] - stats.avg_area;
        var_sum += diff * diff;
    }
    stats.std_dev = sqrtf(var_sum / (float)mesh->triangle_count);

    /* Build histogram */
    if (histogram_bins > 0 && stats.max_area > stats.min_area) {
        stats.histogram = (uint32_t*)calloc(histogram_bins, sizeof(uint32_t));
        stats.bin_width = (stats.max_area - stats.min_area) / (float)histogram_bins;

        for (size_t i = 0; i < mesh->triangle_count; i++) {
            size_t bin = (size_t)((areas[i] - stats.min_area) / stats.bin_width);
            if (bin >= histogram_bins) bin = histogram_bins - 1;
            stats.histogram[bin]++;
        }
    }

    free(areas);
    return stats;
}

void geomlib_area_stats_free(geomlib_area_stats *stats) {
    free(stats->histogram);
}

/* ============================================================================
 * Aspect Ratio Statistics
 * ============================================================================ */

geomlib_aspect_ratio_stats geomlib_compute_aspect_ratio_stats(
    const geomlib_trimesh *mesh,
    size_t histogram_bins
) {
    geomlib_aspect_ratio_stats stats;
    stats.histogram = NULL;
    stats.bin_count = histogram_bins;
    stats.good_count = 0;
    stats.acceptable_count = 0;
    stats.poor_count = 0;

    float *ratios = (float*)malloc(mesh->triangle_count * sizeof(float));
    stats.min_ratio = FLT_MAX;
    stats.max_ratio = 0.0f;
    float sum = 0.0f;

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 v0 = mesh->vertices[mesh->indices[i * 3 + 0]];
        geomlib_vec3 v1 = mesh->vertices[mesh->indices[i * 3 + 1]];
        geomlib_vec3 v2 = mesh->vertices[mesh->indices[i * 3 + 2]];

        float ratio = geomlib_triangle_aspect_ratio(v0, v1, v2);
        ratios[i] = ratio;
        sum += ratio;

        if (ratio < stats.min_ratio) stats.min_ratio = ratio;
        if (ratio > stats.max_ratio) stats.max_ratio = ratio;

        if (ratio < 2.0f) stats.good_count++;
        else if (ratio < 5.0f) stats.acceptable_count++;
        else stats.poor_count++;
    }

    stats.avg_ratio = sum / (float)mesh->triangle_count;

    /* Build histogram */
    if (histogram_bins > 0 && stats.max_ratio > stats.min_ratio) {
        stats.histogram = (uint32_t*)calloc(histogram_bins, sizeof(uint32_t));
        stats.bin_width = (stats.max_ratio - stats.min_ratio) / (float)histogram_bins;

        for (size_t i = 0; i < mesh->triangle_count; i++) {
            size_t bin = (size_t)((ratios[i] - stats.min_ratio) / stats.bin_width);
            if (bin >= histogram_bins) bin = histogram_bins - 1;
            stats.histogram[bin]++;
        }
    }

    free(ratios);
    return stats;
}

void geomlib_aspect_ratio_stats_free(geomlib_aspect_ratio_stats *stats) {
    free(stats->histogram);
}

/* ============================================================================
 * Comprehensive Mesh Statistics
 * ============================================================================ */

geomlib_mesh_stats geomlib_compute_mesh_stats(
    const geomlib_trimesh *mesh,
    size_t histogram_bins
) {
    geomlib_mesh_stats stats;

    stats.vertex_count = mesh->vertex_count;
    stats.triangle_count = mesh->triangle_count;

    /* Compute edge count */
    geomlib_edge_mesh edge_mesh = geomlib_build_edge_mesh(mesh);
    stats.edge_count = edge_mesh.edge_count;

    /* Check if closed */
    stats.is_closed = true;
    for (size_t i = 0; i < edge_mesh.edge_count; i++) {
        if (edge_mesh.edges[i].face1 == UINT32_MAX) {
            stats.is_closed = false;
            break;
        }
    }
    geomlib_edge_mesh_free(&edge_mesh);

    /* Bounding box and centroid */
    stats.bbox = geomlib_trimesh_compute_aabb(mesh);
    stats.centroid = geomlib_trimesh_centroid(mesh);

    /* Sub-statistics */
    stats.edges = geomlib_compute_edge_stats(mesh, histogram_bins);
    stats.angles = geomlib_compute_angle_stats(mesh, histogram_bins);
    stats.areas = geomlib_compute_area_stats(mesh, histogram_bins);
    stats.aspect_ratios = geomlib_compute_aspect_ratio_stats(mesh, histogram_bins);

    stats.surface_area = stats.areas.total_area;
    stats.volume = stats.is_closed ? geomlib_trimesh_volume(mesh) : 0.0f;

    return stats;
}

void geomlib_mesh_stats_free(geomlib_mesh_stats *stats) {
    geomlib_edge_stats_free(&stats->edges);
    geomlib_angle_stats_free(&stats->angles);
    geomlib_area_stats_free(&stats->areas);
    geomlib_aspect_ratio_stats_free(&stats->aspect_ratios);
}

/* ============================================================================
 * Mesh Validation (Basic)
 * ============================================================================ */

geomlib_validation_result geomlib_validate_mesh(
    const geomlib_trimesh *mesh,
    geomlib_validation_level level
) {
    geomlib_validation_result result;
    result.is_valid = true;
    result.level = level;
    result.issues = NULL;
    result.issue_count = 0;
    result.error_count = 0;
    result.warning_count = 0;

    /* TODO: Implement validation */

    return result;
}

void geomlib_validation_result_free(geomlib_validation_result *result) {
    free(result->issues);
}

/* ============================================================================
 * Mesh Comparison (Stub)
 * ============================================================================ */

geomlib_mesh_difference geomlib_compare_meshes(
    const geomlib_trimesh *mesh1,
    const geomlib_trimesh *mesh2
) {
    geomlib_mesh_difference diff;
    memset(&diff, 0, sizeof(diff));
    /* TODO: Implement */
    return diff;
}
