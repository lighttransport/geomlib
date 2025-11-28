/*
 * geomlib - Feature Detection Implementation
 */

#include "geomlib/geomlib_feature.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

static float deg_to_rad(float deg) {
    return deg * (float)M_PI / 180.0f;
}

static float rad_to_deg(float rad) {
    return rad * 180.0f / (float)M_PI;
}

/* ============================================================================
 * Parameters
 * ============================================================================ */

geomlib_feature_params geomlib_feature_params_default(void) {
    geomlib_feature_params params;
    params.sharp_angle_threshold = 30.0f;       /* 30 degrees */
    params.detect_ridges_valleys = true;
    params.curvature_threshold = 0.1f;
    params.detect_corners = true;
    params.min_corner_edges = 3;
    params.mark_boundaries = true;
    params.min_edge_length = 0.0f;              /* No filtering */
    params.filter_small_features = false;
    return params;
}

/* ============================================================================
 * Dihedral Angle Computation
 * ============================================================================ */

float geomlib_compute_dihedral_angle(
    const geomlib_trimesh *mesh,
    uint32_t v0,
    uint32_t v1,
    uint32_t face0,
    uint32_t face1
) {
    if (face1 == UINT32_MAX) {
        return 180.0f;  /* Boundary edge */
    }

    /* Compute face normals */
    geomlib_vec3 t0_v0 = mesh->vertices[mesh->indices[face0 * 3 + 0]];
    geomlib_vec3 t0_v1 = mesh->vertices[mesh->indices[face0 * 3 + 1]];
    geomlib_vec3 t0_v2 = mesh->vertices[mesh->indices[face0 * 3 + 2]];

    geomlib_vec3 e0_0 = geomlib_vec3_sub(t0_v1, t0_v0);
    geomlib_vec3 e0_1 = geomlib_vec3_sub(t0_v2, t0_v0);
    geomlib_vec3 n0 = geomlib_vec3_normalize(geomlib_vec3_cross(e0_0, e0_1));

    geomlib_vec3 t1_v0 = mesh->vertices[mesh->indices[face1 * 3 + 0]];
    geomlib_vec3 t1_v1 = mesh->vertices[mesh->indices[face1 * 3 + 1]];
    geomlib_vec3 t1_v2 = mesh->vertices[mesh->indices[face1 * 3 + 2]];

    geomlib_vec3 e1_0 = geomlib_vec3_sub(t1_v1, t1_v0);
    geomlib_vec3 e1_1 = geomlib_vec3_sub(t1_v2, t1_v0);
    geomlib_vec3 n1 = geomlib_vec3_normalize(geomlib_vec3_cross(e1_0, e1_1));

    /* Dihedral angle is 180 - angle between normals */
    float dot = geomlib_vec3_dot(n0, n1);
    dot = fmaxf(-1.0f, fminf(1.0f, dot));  /* Clamp to [-1, 1] */

    float angle_rad = acosf(dot);
    float angle_deg = rad_to_deg(angle_rad);

    /* Return exterior dihedral angle */
    return 180.0f - angle_deg;
}

/* ============================================================================
 * Sharp Edge Detection
 * ============================================================================ */

bool geomlib_is_edge_sharp(
    const geomlib_trimesh *mesh,
    uint32_t v0,
    uint32_t v1,
    float angle_threshold_deg
) {
    /* Find faces adjacent to edge */
    uint32_t face0 = UINT32_MAX, face1 = UINT32_MAX;
    int found = 0;

    for (size_t i = 0; i < mesh->triangle_count && found < 2; i++) {
        uint32_t i0 = mesh->indices[i * 3 + 0];
        uint32_t i1 = mesh->indices[i * 3 + 1];
        uint32_t i2 = mesh->indices[i * 3 + 2];

        /* Check if edge (v0, v1) is in this triangle */
        bool has_edge = false;
        if ((i0 == v0 && i1 == v1) || (i1 == v0 && i0 == v1) ||
            (i1 == v0 && i2 == v1) || (i2 == v0 && i1 == v1) ||
            (i2 == v0 && i0 == v1) || (i0 == v0 && i2 == v1)) {
            has_edge = true;
        }

        if (has_edge) {
            if (face0 == UINT32_MAX) {
                face0 = (uint32_t)i;
            } else {
                face1 = (uint32_t)i;
            }
            found++;
        }
    }

    if (face0 == UINT32_MAX) {
        return false;  /* Edge not found */
    }

    float dihedral = geomlib_compute_dihedral_angle(mesh, v0, v1, face0, face1);
    return fabsf(dihedral) >= angle_threshold_deg;
}

geomlib_feature_edge* geomlib_detect_sharp_edges(
    const geomlib_trimesh *mesh,
    float angle_threshold_deg,
    size_t *edge_count
) {
    /* Build edge map */
    geomlib_edge_mesh edges = geomlib_build_edge_mesh(mesh);

    /* Allocate space for feature edges */
    geomlib_feature_edge *feature_edges = (geomlib_feature_edge*)malloc(
        edges.edge_count * sizeof(geomlib_feature_edge));
    size_t feature_count = 0;

    /* Check each edge */
    for (size_t i = 0; i < edges.edge_count; i++) {
        uint32_t v0 = edges.edges[i].v0;
        uint32_t v1 = edges.edges[i].v1;
        uint32_t face0 = edges.edges[i].face0;
        uint32_t face1 = edges.edges[i].face1;

        float dihedral = geomlib_compute_dihedral_angle(mesh, v0, v1, face0, face1);

        if (fabsf(dihedral) >= angle_threshold_deg) {
            feature_edges[feature_count].v0 = v0;
            feature_edges[feature_count].v1 = v1;
            feature_edges[feature_count].type = (face1 == UINT32_MAX) ?
                GEOMLIB_FEATURE_BOUNDARY : GEOMLIB_FEATURE_SHARP_EDGE;
            feature_edges[feature_count].sharpness = fabsf(dihedral);
            feature_edges[feature_count].face0 = face0;
            feature_edges[feature_count].face1 = face1;
            feature_count++;
        }
    }

    geomlib_edge_mesh_free(&edges);

    *edge_count = feature_count;
    return feature_edges;
}

/* ============================================================================
 * Boundary Feature Detection
 * ============================================================================ */

geomlib_feature_edge* geomlib_mark_boundary_features(
    const geomlib_trimesh *mesh,
    size_t *edge_count
) {
    geomlib_edge_mesh edges = geomlib_build_edge_mesh(mesh);

    /* Count boundary edges */
    size_t boundary_count = 0;
    for (size_t i = 0; i < edges.edge_count; i++) {
        if (edges.edges[i].face1 == UINT32_MAX) {
            boundary_count++;
        }
    }

    /* Allocate and fill boundary features */
    geomlib_feature_edge *boundary_edges = (geomlib_feature_edge*)malloc(
        boundary_count * sizeof(geomlib_feature_edge));

    size_t idx = 0;
    for (size_t i = 0; i < edges.edge_count; i++) {
        if (edges.edges[i].face1 == UINT32_MAX) {
            boundary_edges[idx].v0 = edges.edges[i].v0;
            boundary_edges[idx].v1 = edges.edges[i].v1;
            boundary_edges[idx].type = GEOMLIB_FEATURE_BOUNDARY;
            boundary_edges[idx].sharpness = 180.0f;
            boundary_edges[idx].face0 = edges.edges[i].face0;
            boundary_edges[idx].face1 = UINT32_MAX;
            idx++;
        }
    }

    geomlib_edge_mesh_free(&edges);

    *edge_count = boundary_count;
    return boundary_edges;
}

/* ============================================================================
 * Ridge and Valley Detection
 * ============================================================================ */

bool geomlib_is_ridge(
    const geomlib_trimesh *mesh,
    const geomlib_mesh_curvature *curvature,
    uint32_t v0,
    uint32_t v1,
    float threshold
) {
    /* Average max principal curvature at edge endpoints */
    float k0 = curvature->curvatures[v0].principal_max;
    float k1 = curvature->curvatures[v1].principal_max;
    float avg_k = (k0 + k1) * 0.5f;

    return avg_k > threshold;
}

bool geomlib_is_valley(
    const geomlib_trimesh *mesh,
    const geomlib_mesh_curvature *curvature,
    uint32_t v0,
    uint32_t v1,
    float threshold
) {
    /* Average min principal curvature at edge endpoints */
    float k0 = curvature->curvatures[v0].principal_min;
    float k1 = curvature->curvatures[v1].principal_min;
    float avg_k = (k0 + k1) * 0.5f;

    return avg_k < -threshold;
}

geomlib_feature_edge* geomlib_detect_ridges_valleys(
    const geomlib_trimesh *mesh,
    const geomlib_mesh_curvature *curvature,
    float threshold,
    size_t *edge_count
) {
    geomlib_edge_mesh edges = geomlib_build_edge_mesh(mesh);

    geomlib_feature_edge *feature_edges = (geomlib_feature_edge*)malloc(
        edges.edge_count * sizeof(geomlib_feature_edge));
    size_t feature_count = 0;

    for (size_t i = 0; i < edges.edge_count; i++) {
        uint32_t v0 = edges.edges[i].v0;
        uint32_t v1 = edges.edges[i].v1;

        if (geomlib_is_ridge(mesh, curvature, v0, v1, threshold)) {
            feature_edges[feature_count].v0 = v0;
            feature_edges[feature_count].v1 = v1;
            feature_edges[feature_count].type = GEOMLIB_FEATURE_RIDGE;
            feature_edges[feature_count].sharpness =
                (curvature->curvatures[v0].principal_max + curvature->curvatures[v1].principal_max) * 0.5f;
            feature_edges[feature_count].face0 = edges.edges[i].face0;
            feature_edges[feature_count].face1 = edges.edges[i].face1;
            feature_count++;
        } else if (geomlib_is_valley(mesh, curvature, v0, v1, threshold)) {
            feature_edges[feature_count].v0 = v0;
            feature_edges[feature_count].v1 = v1;
            feature_edges[feature_count].type = GEOMLIB_FEATURE_VALLEY;
            feature_edges[feature_count].sharpness =
                fabsf((curvature->curvatures[v0].principal_min + curvature->curvatures[v1].principal_min) * 0.5f);
            feature_edges[feature_count].face0 = edges.edges[i].face0;
            feature_edges[feature_count].face1 = edges.edges[i].face1;
            feature_count++;
        }
    }

    geomlib_edge_mesh_free(&edges);

    *edge_count = feature_count;
    return feature_edges;
}

/* ============================================================================
 * Corner Detection
 * ============================================================================ */

uint32_t geomlib_count_feature_edges_at_vertex(
    const geomlib_feature_edge *edges,
    size_t edge_count,
    uint32_t vertex_index
) {
    uint32_t count = 0;
    for (size_t i = 0; i < edge_count; i++) {
        if (edges[i].v0 == vertex_index || edges[i].v1 == vertex_index) {
            count++;
        }
    }
    return count;
}

void geomlib_get_feature_edges_at_vertex(
    const geomlib_feature_edge *edges,
    size_t edge_count,
    uint32_t vertex_index,
    geomlib_feature_edge **vertex_edges,
    size_t *vertex_edge_count
) {
    uint32_t count = geomlib_count_feature_edges_at_vertex(edges, edge_count, vertex_index);

    *vertex_edges = (geomlib_feature_edge*)malloc(count * sizeof(geomlib_feature_edge));
    *vertex_edge_count = count;

    size_t idx = 0;
    for (size_t i = 0; i < edge_count; i++) {
        if (edges[i].v0 == vertex_index || edges[i].v1 == vertex_index) {
            (*vertex_edges)[idx++] = edges[i];
        }
    }
}

bool geomlib_is_corner(
    const geomlib_feature_edge *feature_edges,
    size_t edge_count,
    uint32_t vertex_index,
    uint32_t min_edges
) {
    uint32_t count = geomlib_count_feature_edges_at_vertex(
        feature_edges, edge_count, vertex_index);
    return count >= min_edges;
}

geomlib_feature_vertex* geomlib_detect_corners(
    const geomlib_trimesh *mesh,
    const geomlib_feature_edge *feature_edges,
    size_t edge_count,
    uint32_t min_edges,
    size_t *corner_count
) {
    /* Find unique vertices in feature edges */
    bool *is_feature_vertex = (bool*)calloc(mesh->vertex_count, sizeof(bool));
    uint32_t *edge_counts = (uint32_t*)calloc(mesh->vertex_count, sizeof(uint32_t));

    for (size_t i = 0; i < edge_count; i++) {
        is_feature_vertex[feature_edges[i].v0] = true;
        is_feature_vertex[feature_edges[i].v1] = true;
        edge_counts[feature_edges[i].v0]++;
        edge_counts[feature_edges[i].v1]++;
    }

    /* Count corners */
    size_t count = 0;
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        if (is_feature_vertex[i] && edge_counts[i] >= min_edges) {
            count++;
        }
    }

    /* Allocate and fill corners */
    geomlib_feature_vertex *corners = (geomlib_feature_vertex*)malloc(
        count * sizeof(geomlib_feature_vertex));

    size_t idx = 0;
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        if (is_feature_vertex[i] && edge_counts[i] >= min_edges) {
            corners[idx].vertex_index = (uint32_t)i;
            corners[idx].type = GEOMLIB_FEATURE_CORNER;
            corners[idx].edge_count = edge_counts[i];

            /* Find max sharpness */
            float max_sharp = 0.0f;
            for (size_t j = 0; j < edge_count; j++) {
                if (feature_edges[j].v0 == i || feature_edges[j].v1 == i) {
                    if (feature_edges[j].sharpness > max_sharp) {
                        max_sharp = feature_edges[j].sharpness;
                    }
                }
            }
            corners[idx].sharpness = max_sharp;
            idx++;
        }
    }

    free(is_feature_vertex);
    free(edge_counts);

    *corner_count = count;
    return corners;
}

/* ============================================================================
 * Feature Curve Extraction
 * ============================================================================ */

int geomlib_trace_feature_curve(
    const geomlib_feature_edge *edges,
    size_t edge_count,
    uint32_t start_vertex,
    bool *visited_edges,
    uint32_t *curve,
    size_t *curve_length
) {
    size_t length = 0;
    uint32_t current = start_vertex;
    curve[length++] = current;

    bool found = true;
    while (found) {
        found = false;

        /* Find unvisited edge connected to current vertex */
        for (size_t i = 0; i < edge_count; i++) {
            if (visited_edges[i]) continue;

            uint32_t next = UINT32_MAX;
            if (edges[i].v0 == current) {
                next = edges[i].v1;
            } else if (edges[i].v1 == current) {
                next = edges[i].v0;
            }

            if (next != UINT32_MAX) {
                visited_edges[i] = true;
                curve[length++] = next;
                current = next;
                found = true;
                break;
            }
        }
    }

    *curve_length = length;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_extract_feature_curves(
    const geomlib_feature_edge *edges,
    size_t edge_count,
    uint32_t ***curves,
    size_t **curve_lengths,
    size_t *curve_count
) {
    if (edge_count == 0) {
        *curves = NULL;
        *curve_lengths = NULL;
        *curve_count = 0;
        return GEOMLIB_MESH_SUCCESS;
    }

    bool *visited = (bool*)calloc(edge_count, sizeof(bool));

    /* Temporary storage */
    uint32_t **temp_curves = (uint32_t**)malloc(edge_count * sizeof(uint32_t*));
    size_t *temp_lengths = (size_t*)malloc(edge_count * sizeof(size_t));
    size_t num_curves = 0;

    /* Start from each unvisited edge */
    for (size_t i = 0; i < edge_count; i++) {
        if (visited[i]) continue;

        uint32_t *curve = (uint32_t*)malloc(edge_count * 2 * sizeof(uint32_t));
        size_t length;

        visited[i] = true;
        geomlib_trace_feature_curve(edges, edge_count, edges[i].v0,
                                      visited, curve, &length);

        temp_curves[num_curves] = curve;
        temp_lengths[num_curves] = length;
        num_curves++;
    }

    free(visited);

    *curves = temp_curves;
    *curve_lengths = temp_lengths;
    *curve_count = num_curves;

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Main Feature Detection
 * ============================================================================ */

geomlib_features geomlib_detect_features(
    const geomlib_trimesh *mesh,
    const geomlib_feature_params *params
) {
    geomlib_features features;
    memset(&features, 0, sizeof(geomlib_features));

    /* Detect sharp edges */
    size_t sharp_count = 0;
    geomlib_feature_edge *sharp_edges = geomlib_detect_sharp_edges(
        mesh, params->sharp_angle_threshold, &sharp_count);

    /* Detect boundaries if requested */
    size_t boundary_count = 0;
    geomlib_feature_edge *boundary_edges = NULL;
    if (params->mark_boundaries) {
        boundary_edges = geomlib_mark_boundary_features(mesh, &boundary_count);
    }

    /* Combine all edges */
    size_t total_edges = sharp_count + boundary_count;
    features.edges = (geomlib_feature_edge*)malloc(
        total_edges * sizeof(geomlib_feature_edge));
    features.edge_count = total_edges;

    memcpy(features.edges, sharp_edges, sharp_count * sizeof(geomlib_feature_edge));
    if (boundary_count > 0) {
        memcpy(features.edges + sharp_count, boundary_edges,
               boundary_count * sizeof(geomlib_feature_edge));
    }

    free(sharp_edges);
    free(boundary_edges);

    /* Detect corners if requested */
    if (params->detect_corners) {
        features.vertices = geomlib_detect_corners(
            mesh, features.edges, features.edge_count,
            params->min_corner_edges, &features.vertex_count);
    }

    /* Extract feature curves */
    geomlib_extract_feature_curves(
        features.edges, features.edge_count,
        &features.curves, &features.curve_lengths, &features.curve_count);

    return features;
}

void geomlib_features_free(geomlib_features *features) {
    free(features->edges);
    free(features->vertices);

    for (size_t i = 0; i < features->curve_count; i++) {
        free(features->curves[i]);
    }
    free(features->curves);
    free(features->curve_lengths);

    memset(features, 0, sizeof(geomlib_features));
}

/* ============================================================================
 * Planar Region Segmentation
 * ============================================================================ */

geomlib_regions geomlib_segment_planar_regions(
    const geomlib_trimesh *mesh,
    float angle_threshold_deg,
    float distance_threshold
) {
    geomlib_regions regions;
    memset(&regions, 0, sizeof(geomlib_regions));

    /* Simple region growing */
    regions.triangle_labels = (uint32_t*)malloc(mesh->triangle_count * sizeof(uint32_t));
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        regions.triangle_labels[i] = UINT32_MAX;  /* Unlabeled */
    }

    /* Allocate temp regions */
    geomlib_planar_region *temp_regions = (geomlib_planar_region*)malloc(
        mesh->triangle_count * sizeof(geomlib_planar_region));
    size_t region_count = 0;

    /* TODO: Implement proper region growing algorithm */
    /* For now, just return empty regions */

    regions.regions = temp_regions;
    regions.region_count = region_count;

    return regions;
}

void geomlib_regions_free(geomlib_regions *regions) {
    for (size_t i = 0; i < regions->region_count; i++) {
        free(regions->regions[i].triangles);
    }
    free(regions->regions);
    free(regions->triangle_labels);
    memset(regions, 0, sizeof(geomlib_regions));
}

/* ============================================================================
 * Feature Statistics
 * ============================================================================ */

geomlib_feature_stats geomlib_compute_feature_stats(
    const geomlib_features *features
) {
    geomlib_feature_stats stats;
    memset(&stats, 0, sizeof(geomlib_feature_stats));

    stats.corner_count = features->vertex_count;

    float sum_sharpness = 0.0f;
    stats.max_sharpness = -FLT_MAX;
    stats.min_sharpness = FLT_MAX;

    for (size_t i = 0; i < features->edge_count; i++) {
        switch (features->edges[i].type) {
            case GEOMLIB_FEATURE_SHARP_EDGE:
                stats.sharp_edge_count++;
                break;
            case GEOMLIB_FEATURE_RIDGE:
                stats.ridge_count++;
                break;
            case GEOMLIB_FEATURE_VALLEY:
                stats.valley_count++;
                break;
            case GEOMLIB_FEATURE_BOUNDARY:
                stats.boundary_edge_count++;
                break;
            default:
                break;
        }

        sum_sharpness += features->edges[i].sharpness;
        if (features->edges[i].sharpness > stats.max_sharpness) {
            stats.max_sharpness = features->edges[i].sharpness;
        }
        if (features->edges[i].sharpness < stats.min_sharpness) {
            stats.min_sharpness = features->edges[i].sharpness;
        }
    }

    if (features->edge_count > 0) {
        stats.avg_sharpness = sum_sharpness / (float)features->edge_count;
    }

    return stats;
}
