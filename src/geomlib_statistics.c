/*
 * geomlib - Mesh Statistics Implementation (Basic)
 */

#include "geomlib/geomlib_statistics.h"
#include "geomlib/geomlib_topology.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
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
    /* Sort lengths for median calculation */
    for (size_t i = 0; i < edge_mesh.edge_count - 1; i++) {
        for (size_t j = i + 1; j < edge_mesh.edge_count; j++) {
            if (lengths[i] > lengths[j]) {
                float temp = lengths[i];
                lengths[i] = lengths[j];
                lengths[j] = temp;
            }
        }
    }

    /* Compute median */
    if (edge_mesh.edge_count % 2 == 0) {
        size_t mid = edge_mesh.edge_count / 2;
        stats.median_length = (lengths[mid - 1] + lengths[mid]) / 2.0f;
    } else {
        stats.median_length = lengths[edge_mesh.edge_count / 2];
    }

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

    if (!mesh) {
        return result;
    }

    /* Allocate space for issues (worst case: all triangles have issues) */
    size_t max_issues = mesh->triangle_count * 3;
    geomlib_validation_issue *issues = (geomlib_validation_issue*)malloc(
        max_issues * sizeof(geomlib_validation_issue));
    if (!issues) {
        return result;
    }

    size_t issue_count = 0;

    /* Check for degenerate triangles */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t i0 = mesh->indices[i * 3 + 0];
        uint32_t i1 = mesh->indices[i * 3 + 1];
        uint32_t i2 = mesh->indices[i * 3 + 2];

        /* Check for invalid indices */
        if (i0 >= mesh->vertex_count || i1 >= mesh->vertex_count || i2 >= mesh->vertex_count) {
            issues[issue_count].type = GEOMLIB_ISSUE_INVALID_INDEX;
            issues[issue_count].element_index = (uint32_t)i;
            snprintf(issues[issue_count].message, 256,
                    "Triangle %zu has out-of-bounds vertex index", i);
            issue_count++;
            result.error_count++;
            result.is_valid = false;
            continue;
        }

        /* Check for duplicate vertices in triangle */
        if (i0 == i1 || i1 == i2 || i2 == i0) {
            issues[issue_count].type = GEOMLIB_ISSUE_DUPLICATE_VERTEX;
            issues[issue_count].element_index = (uint32_t)i;
            snprintf(issues[issue_count].message, 256,
                    "Triangle %zu has duplicate vertices", i);
            issue_count++;
            result.error_count++;
            result.is_valid = false;
            continue;
        }

        /* Check for zero-area triangles */
        geomlib_vec3 v0 = mesh->vertices[i0];
        geomlib_vec3 v1 = mesh->vertices[i1];
        geomlib_vec3 v2 = mesh->vertices[i2];

        geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
        geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
        geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);
        float area = geomlib_vec3_length(cross) * 0.5f;

        if (area < 1e-10f) {
            issues[issue_count].type = GEOMLIB_ISSUE_ZERO_AREA_TRIANGLE;
            issues[issue_count].element_index = (uint32_t)i;
            snprintf(issues[issue_count].message, 256,
                    "Triangle %zu is degenerate (near-zero area)", i);
            issue_count++;
            result.warning_count++;
            if (level >= GEOMLIB_VALIDATE_STRICT) {
                result.is_valid = false;
            }
        }
    }

    /* Check for NaN or infinite vertices */
    if (level >= GEOMLIB_VALIDATE_STRICT) {
        for (size_t i = 0; i < mesh->vertex_count; i++) {
            geomlib_vec3 v = mesh->vertices[i];
            if (isnan(v.x) || isnan(v.y) || isnan(v.z) ||
                isinf(v.x) || isinf(v.y) || isinf(v.z)) {
                issues[issue_count].type = GEOMLIB_ISSUE_INVALID_INDEX;
                issues[issue_count].element_index = (uint32_t)i;
                snprintf(issues[issue_count].message, 256,
                        "Vertex %zu has NaN or infinite coordinates", i);
                issue_count++;
                result.error_count++;
                result.is_valid = false;
            }
        }
    }

    /* Store results */
    if (issue_count > 0) {
        result.issues = (geomlib_validation_issue*)realloc(issues,
            issue_count * sizeof(geomlib_validation_issue));
        result.issue_count = issue_count;
    } else {
        free(issues);
    }

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

    if (!mesh1 || !mesh2) {
        return diff;
    }

    /* Vertex and triangle count differences */
    diff.vertex_count_diff = (mesh1->vertex_count > mesh2->vertex_count) ?
        (mesh1->vertex_count - mesh2->vertex_count) :
        (mesh2->vertex_count - mesh1->vertex_count);

    diff.triangle_count_diff = (mesh1->triangle_count > mesh2->triangle_count) ?
        (mesh1->triangle_count - mesh2->triangle_count) :
        (mesh2->triangle_count - mesh1->triangle_count);

    /* Compute Hausdorff distance (one direction: mesh1 to mesh2) */
    float max_dist_1to2 = 0.0f;
    float sum_dist_1to2 = 0.0f;
    float sum_sq_dist_1to2 = 0.0f;

    for (size_t i = 0; i < mesh1->vertex_count; i++) {
        geomlib_vec3 v1 = mesh1->vertices[i];
        float min_dist = FLT_MAX;

        /* Find closest point on mesh2 */
        for (size_t j = 0; j < mesh2->triangle_count; j++) {
            geomlib_vec3 t0 = mesh2->vertices[mesh2->indices[j * 3 + 0]];
            geomlib_vec3 t1 = mesh2->vertices[mesh2->indices[j * 3 + 1]];
            geomlib_vec3 t2 = mesh2->vertices[mesh2->indices[j * 3 + 2]];

            /* Simple vertex-to-vertex distance for now */
            /* TODO: Proper point-to-triangle distance */
            float d0 = geomlib_vec3_length(geomlib_vec3_sub(v1, t0));
            float d1 = geomlib_vec3_length(geomlib_vec3_sub(v1, t1));
            float d2 = geomlib_vec3_length(geomlib_vec3_sub(v1, t2));

            float d_min = fminf(fminf(d0, d1), d2);
            if (d_min < min_dist) {
                min_dist = d_min;
            }
        }

        if (min_dist > max_dist_1to2) {
            max_dist_1to2 = min_dist;
        }
        sum_dist_1to2 += min_dist;
        sum_sq_dist_1to2 += min_dist * min_dist;
    }

    /* Compute Hausdorff distance (other direction: mesh2 to mesh1) */
    float max_dist_2to1 = 0.0f;
    float sum_dist_2to1 = 0.0f;
    float sum_sq_dist_2to1 = 0.0f;

    for (size_t i = 0; i < mesh2->vertex_count; i++) {
        geomlib_vec3 v2 = mesh2->vertices[i];
        float min_dist = FLT_MAX;

        /* Find closest point on mesh1 */
        for (size_t j = 0; j < mesh1->triangle_count; j++) {
            geomlib_vec3 t0 = mesh1->vertices[mesh1->indices[j * 3 + 0]];
            geomlib_vec3 t1 = mesh1->vertices[mesh1->indices[j * 3 + 1]];
            geomlib_vec3 t2 = mesh1->vertices[mesh1->indices[j * 3 + 2]];

            float d0 = geomlib_vec3_length(geomlib_vec3_sub(v2, t0));
            float d1 = geomlib_vec3_length(geomlib_vec3_sub(v2, t1));
            float d2 = geomlib_vec3_length(geomlib_vec3_sub(v2, t2));

            float d_min = fminf(fminf(d0, d1), d2);
            if (d_min < min_dist) {
                min_dist = d_min;
            }
        }

        if (min_dist > max_dist_2to1) {
            max_dist_2to1 = min_dist;
        }
        sum_dist_2to1 += min_dist;
        sum_sq_dist_2to1 += min_dist * min_dist;
    }

    /* Hausdorff distance is the maximum of both directions */
    diff.hausdorff_distance = fmaxf(max_dist_1to2, max_dist_2to1);

    /* Average distance */
    size_t total_vertices = mesh1->vertex_count + mesh2->vertex_count;
    diff.avg_distance = (sum_dist_1to2 + sum_dist_2to1) / (float)total_vertices;

    /* RMS distance */
    diff.rms_distance = sqrtf((sum_sq_dist_1to2 + sum_sq_dist_2to1) / (float)total_vertices);

    /* Volume difference (if both meshes are closed) */
    float vol1 = geomlib_trimesh_volume(mesh1);
    float vol2 = geomlib_trimesh_volume(mesh2);
    diff.volume_diff = fabsf(vol1 - vol2);
    if (fabsf(vol1) > 1e-10f) {
        diff.volume_diff_percent = (diff.volume_diff / fabsf(vol1)) * 100.0f;
    } else {
        diff.volume_diff_percent = 0.0f;
    }

    return diff;
}
