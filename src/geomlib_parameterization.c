/*
 * geomlib - Mesh Parameterization Implementation
 * Based on LSCM (Lévy et al. 2002) and harmonic mapping
 */

#include "geomlib/geomlib_parameterization.h"
#include "geomlib/geomlib_healing.h"
#include "geomlib/geomlib_solver.h"
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Helper: Find boundary loop using half-edge structure */
int geomlib_find_boundary_loop(
    const geomlib_trimesh *mesh,
    uint32_t **boundary_vertices,
    size_t *boundary_count
) {
    if (!mesh || !boundary_vertices || !boundary_count) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Build half-edge structure */
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Find boundary half-edges (those without twins) */
    uint32_t *boundary = malloc(mesh->vertex_count * sizeof(uint32_t));
    size_t count = 0;

    /* Find a boundary half-edge */
    uint32_t start_he = UINT32_MAX;
    for (size_t i = 0; i < hemesh.halfedge_count; i++) {
        if (hemesh.halfedges[i].twin == UINT32_MAX) {
            start_he = i;
            break;
        }
    }

    if (start_he == UINT32_MAX) {
        /* No boundary found - closed mesh */
        free(boundary);
        geomlib_halfedge_free(&hemesh);
        *boundary_vertices = NULL;
        *boundary_count = 0;
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Traverse boundary loop */
    uint32_t current_he = start_he;
    do {
        uint32_t vertex = hemesh.halfedges[current_he].vertex;
        boundary[count++] = vertex;

        /* Move to next boundary half-edge */
        uint32_t next_he = hemesh.halfedges[current_he].next;

        /* Find the boundary edge that starts from the end of current edge */
        uint32_t next_vertex = (next_he != UINT32_MAX) ?
                               hemesh.halfedges[next_he].vertex : UINT32_MAX;

        if (next_vertex == UINT32_MAX) break;

        /* Find boundary edge from next_vertex */
        current_he = UINT32_MAX;
        for (size_t i = 0; i < hemesh.halfedge_count; i++) {
            if (hemesh.halfedges[i].vertex == next_vertex &&
                hemesh.halfedges[i].twin == UINT32_MAX) {
                current_he = i;
                break;
            }
        }

        if (current_he == UINT32_MAX) break;
    } while (current_he != start_he && count < mesh->vertex_count);

    geomlib_halfedge_free(&hemesh);

    *boundary_vertices = boundary;
    *boundary_count = count;

    return GEOMLIB_MESH_SUCCESS;
}

/* Pin boundary vertices to unit circle */
void geomlib_pin_boundary_circle(
    const uint32_t *boundary_vertices,
    size_t boundary_count,
    geomlib_vec2 *uvs
) {
    for (size_t i = 0; i < boundary_count; i++) {
        float angle = 2.0f * M_PI * (float)i / boundary_count;
        uint32_t v = boundary_vertices[i];
        uvs[v].x = cosf(angle) * 0.5f + 0.5f;  /* Map to [0,1] */
        uvs[v].y = sinf(angle) * 0.5f + 0.5f;
    }
}

/* Pin boundary vertices to unit square */
void geomlib_pin_boundary_square(
    const uint32_t *boundary_vertices,
    size_t boundary_count,
    geomlib_vec2 *uvs
) {
    /* Distribute vertices along square perimeter */
    for (size_t i = 0; i < boundary_count; i++) {
        float t = (float)i / boundary_count;  /* Parameter [0,1] */
        uint32_t v = boundary_vertices[i];

        if (t < 0.25f) {
            /* Bottom edge */
            uvs[v].x = t * 4.0f;
            uvs[v].y = 0.0f;
        } else if (t < 0.5f) {
            /* Right edge */
            uvs[v].x = 1.0f;
            uvs[v].y = (t - 0.25f) * 4.0f;
        } else if (t < 0.75f) {
            /* Top edge */
            uvs[v].x = 1.0f - (t - 0.5f) * 4.0f;
            uvs[v].y = 1.0f;
        } else {
            /* Left edge */
            uvs[v].x = 0.0f;
            uvs[v].y = 1.0f - (t - 0.75f) * 4.0f;
        }
    }
}

/* Harmonic/Tutte parameterization - simple iterative version */
int geomlib_parameterize_harmonic(
    const geomlib_trimesh *mesh,
    geomlib_boundary_shape boundary,
    geomlib_uv_coords *uvs
) {
    if (!mesh || !uvs) return GEOMLIB_MESH_ERROR_INVALID;

    /* Find boundary */
    uint32_t *boundary_verts = NULL;
    size_t boundary_count = 0;
    int result = geomlib_find_boundary_loop(mesh, &boundary_verts, &boundary_count);
    if (result != GEOMLIB_MESH_SUCCESS || boundary_count == 0) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Initialize UVs */
    memset(uvs->uvs, 0, mesh->vertex_count * sizeof(geomlib_vec2));

    /* Pin boundary */
    if (boundary == GEOMLIB_BOUNDARY_CIRCLE) {
        geomlib_pin_boundary_circle(boundary_verts, boundary_count, uvs->uvs);
    } else {
        geomlib_pin_boundary_square(boundary_verts, boundary_count, uvs->uvs);
    }

    /* Mark boundary vertices */
    bool *is_boundary = calloc(mesh->vertex_count, sizeof(bool));
    for (size_t i = 0; i < boundary_count; i++) {
        is_boundary[boundary_verts[i]] = true;
    }

    /* Build half-edge structure for neighbor queries */
    geomlib_halfedge_mesh hemesh;
    result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) {
        free(boundary_verts);
        free(is_boundary);
        return result;
    }

    /* Solve using simple Jacobi iteration */
    /* Each interior vertex = average of its neighbors */
    geomlib_vec2 *uvs_new = malloc(mesh->vertex_count * sizeof(geomlib_vec2));

    const int max_iterations = 1000;
    const float tolerance = 1e-6f;

    for (int iter = 0; iter < max_iterations; iter++) {
        float max_change = 0.0f;

        /* Update each interior vertex */
        for (size_t v = 0; v < mesh->vertex_count; v++) {
            if (is_boundary[v]) {
                uvs_new[v] = uvs->uvs[v];
                continue;
            }

            /* Compute average of neighbors */
            geomlib_vec2 avg = {0.0f, 0.0f};
            int neighbor_count = 0;

            /* Find one-ring neighbors via half-edge */
            for (size_t he_idx = 0; he_idx < hemesh.halfedge_count; he_idx++) {
                if (hemesh.halfedges[he_idx].vertex == v) {
                    uint32_t next = hemesh.halfedges[he_idx].next;
                    if (next != UINT32_MAX) {
                        uint32_t neighbor = hemesh.halfedges[next].vertex;
                        avg.x += uvs->uvs[neighbor].x;
                        avg.y += uvs->uvs[neighbor].y;
                        neighbor_count++;
                    }
                }
            }

            if (neighbor_count > 0) {
                avg.x /= neighbor_count;
                avg.y /= neighbor_count;
            }

            uvs_new[v] = avg;

            /* Compute change */
            float dx = avg.x - uvs->uvs[v].x;
            float dy = avg.y - uvs->uvs[v].y;
            float change = sqrtf(dx*dx + dy*dy);
            if (change > max_change) max_change = change;
        }

        /* Copy new values */
        memcpy(uvs->uvs, uvs_new, mesh->vertex_count * sizeof(geomlib_vec2));

        /* Check convergence */
        if (max_change < tolerance) {
            break;
        }
    }

    /* Cleanup */
    free(uvs_new);
    free(is_boundary);
    free(boundary_verts);
    geomlib_halfedge_free(&hemesh);

    return GEOMLIB_MESH_SUCCESS;
}

/* Helper: Compute LSCM conformal weights for a triangle */
static void geomlib_compute_lscm_weights(
    geomlib_vec3 p0, geomlib_vec3 p1, geomlib_vec3 p2,
    float *w0, float *w1, float *w2
) {
    /* Compute edges */
    geomlib_vec3 e0 = geomlib_vec3_sub(p1, p0);
    geomlib_vec3 e1 = geomlib_vec3_sub(p2, p0);
    geomlib_vec3 e2 = geomlib_vec3_sub(p2, p1);

    /* Compute triangle area */
    geomlib_vec3 n = geomlib_vec3_cross(e0, e1);
    float area = geomlib_vec3_length(n) * 0.5f;

    if (area < 1e-10f) {
        *w0 = *w1 = *w2 = 0.0f;
        return;
    }

    /* Compute cotangent weights (used in conformal parameterization) */
    float len_e0 = geomlib_vec3_length(e0);
    float len_e1 = geomlib_vec3_length(e1);
    float len_e2 = geomlib_vec3_length(e2);

    if (len_e0 < 1e-10f || len_e1 < 1e-10f || len_e2 < 1e-10f) {
        *w0 = *w1 = *w2 = 0.0f;
        return;
    }

    /* Angles using law of cosines */
    float cos_alpha = geomlib_vec3_dot(e0, e1) / (len_e0 * len_e1);
    float cos_beta = geomlib_vec3_dot(geomlib_vec3_negate(e0), e2) / (len_e0 * len_e2);
    float cos_gamma = geomlib_vec3_dot(geomlib_vec3_negate(e1), geomlib_vec3_negate(e2)) / (len_e1 * len_e2);

    /* Clamp to prevent numerical issues */
    cos_alpha = fmaxf(-1.0f, fminf(1.0f, cos_alpha));
    cos_beta = fmaxf(-1.0f, fminf(1.0f, cos_beta));
    cos_gamma = fmaxf(-1.0f, fminf(1.0f, cos_gamma));

    float sin_alpha = sqrtf(1.0f - cos_alpha * cos_alpha);
    float sin_beta = sqrtf(1.0f - cos_beta * cos_beta);
    float sin_gamma = sqrtf(1.0f - cos_gamma * cos_gamma);

    /* Cotangent weights for conformal mapping */
    *w0 = (sin_beta > 1e-6f) ? (cos_beta / sin_beta) : 0.0f;
    *w1 = (sin_gamma > 1e-6f) ? (cos_gamma / sin_gamma) : 0.0f;
    *w2 = (sin_alpha > 1e-6f) ? (cos_alpha / sin_alpha) : 0.0f;
}

/* LSCM parameterization (Lévy et al. 2002) */
int geomlib_parameterize_lscm(
    const geomlib_trimesh *mesh,
    geomlib_boundary_shape boundary,
    geomlib_uv_coords *uvs
) {
    if (!mesh || !uvs) return GEOMLIB_MESH_ERROR_INVALID;

    /* Find boundary */
    uint32_t *boundary_verts = NULL;
    size_t boundary_count = 0;
    int result = geomlib_find_boundary_loop(mesh, &boundary_verts, &boundary_count);
    if (result != GEOMLIB_MESH_SUCCESS || boundary_count == 0) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Initialize UVs */
    memset(uvs->uvs, 0, mesh->vertex_count * sizeof(geomlib_vec2));

    /* Pin boundary */
    if (boundary == GEOMLIB_BOUNDARY_CIRCLE) {
        geomlib_pin_boundary_circle(boundary_verts, boundary_count, uvs->uvs);
    } else {
        geomlib_pin_boundary_square(boundary_verts, boundary_count, uvs->uvs);
    }

    /* Mark boundary vertices */
    bool *is_boundary = calloc(mesh->vertex_count, sizeof(bool));
    for (size_t i = 0; i < boundary_count; i++) {
        is_boundary[boundary_verts[i]] = true;
    }

    /* Build half-edge structure */
    geomlib_halfedge_mesh hemesh;
    result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) {
        free(boundary_verts);
        free(is_boundary);
        return result;
    }

    /* LSCM iterative solver using conformal weights */
    geomlib_vec2 *uvs_new = malloc(mesh->vertex_count * sizeof(geomlib_vec2));

    const int max_iterations = 1000;
    const float tolerance = 1e-6f;

    for (int iter = 0; iter < max_iterations; iter++) {
        float max_change = 0.0f;

        /* Update each interior vertex using conformal Laplacian */
        for (size_t v = 0; v < mesh->vertex_count; v++) {
            if (is_boundary[v]) {
                uvs_new[v] = uvs->uvs[v];
                continue;
            }

            /* Compute weighted average using cotangent weights */
            geomlib_vec2 weighted_sum = {0.0f, 0.0f};
            float weight_sum = 0.0f;

            /* Iterate over triangles incident to this vertex */
            for (size_t t = 0; t < mesh->triangle_count; t++) {
                uint32_t i0 = mesh->indices[t * 3 + 0];
                uint32_t i1 = mesh->indices[t * 3 + 1];
                uint32_t i2 = mesh->indices[t * 3 + 2];

                /* Check if vertex v is in this triangle */
                int vertex_pos = -1;
                if (i0 == v) vertex_pos = 0;
                else if (i1 == v) vertex_pos = 1;
                else if (i2 == v) vertex_pos = 2;

                if (vertex_pos < 0) continue;

                /* Get triangle vertices */
                geomlib_vec3 p0 = mesh->vertices[i0];
                geomlib_vec3 p1 = mesh->vertices[i1];
                geomlib_vec3 p2 = mesh->vertices[i2];

                /* Compute conformal weights */
                float w0, w1, w2;
                geomlib_compute_lscm_weights(p0, p1, p2, &w0, &w1, &w2);

                /* Add weighted contributions from neighbors */
                if (vertex_pos == 0) {
                    /* Vertex is p0, neighbors are p1 and p2 */
                    weighted_sum.x += w1 * uvs->uvs[i1].x + w2 * uvs->uvs[i2].x;
                    weighted_sum.y += w1 * uvs->uvs[i1].y + w2 * uvs->uvs[i2].y;
                    weight_sum += w1 + w2;
                } else if (vertex_pos == 1) {
                    /* Vertex is p1, neighbors are p0 and p2 */
                    weighted_sum.x += w2 * uvs->uvs[i0].x + w0 * uvs->uvs[i2].x;
                    weighted_sum.y += w2 * uvs->uvs[i0].y + w0 * uvs->uvs[i2].y;
                    weight_sum += w2 + w0;
                } else {
                    /* Vertex is p2, neighbors are p0 and p1 */
                    weighted_sum.x += w0 * uvs->uvs[i0].x + w1 * uvs->uvs[i1].x;
                    weighted_sum.y += w0 * uvs->uvs[i0].y + w1 * uvs->uvs[i1].y;
                    weight_sum += w0 + w1;
                }
            }

            /* Compute new position */
            if (weight_sum > 1e-10f) {
                uvs_new[v].x = weighted_sum.x / weight_sum;
                uvs_new[v].y = weighted_sum.y / weight_sum;
            } else {
                uvs_new[v] = uvs->uvs[v];
            }

            /* Compute change */
            float dx = uvs_new[v].x - uvs->uvs[v].x;
            float dy = uvs_new[v].y - uvs->uvs[v].y;
            float change = sqrtf(dx*dx + dy*dy);
            if (change > max_change) max_change = change;
        }

        /* Copy new values */
        memcpy(uvs->uvs, uvs_new, mesh->vertex_count * sizeof(geomlib_vec2));

        /* Check convergence */
        if (max_change < tolerance) {
            break;
        }
    }

    /* Cleanup */
    free(uvs_new);
    free(is_boundary);
    free(boundary_verts);
    geomlib_halfedge_free(&hemesh);

    return GEOMLIB_MESH_SUCCESS;
}

/* Default parameters */
geomlib_param_params geomlib_param_params_default(geomlib_param_method method) {
    geomlib_param_params params;
    params.method = method;
    params.boundary = GEOMLIB_BOUNDARY_CIRCLE;
    params.pin_vertex_0 = UINT32_MAX;
    params.pin_vertex_1 = UINT32_MAX;
    params.max_iterations = 1000;
    params.tolerance = 1e-6f;
    params.normalize_uvs = true;
    return params;
}

/* Main parameterization function */
int geomlib_parameterize(
    const geomlib_trimesh *mesh,
    const geomlib_param_params *params,
    geomlib_uv_coords *uvs
) {
    if (!mesh || !params || !uvs) return GEOMLIB_MESH_ERROR_INVALID;

    int result;
    switch (params->method) {
        case GEOMLIB_PARAM_LSCM:
            result = geomlib_parameterize_lscm(mesh, params->boundary, uvs);
            break;
        case GEOMLIB_PARAM_HARMONIC:
            result = geomlib_parameterize_harmonic(mesh, params->boundary, uvs);
            break;
        default:
            return GEOMLIB_MESH_ERROR_INVALID;
    }

    if (result == GEOMLIB_MESH_SUCCESS && params->normalize_uvs) {
        geomlib_normalize_uvs(uvs);
    }

    return result;
}

/* Normalize UVs to [0,1] range */
void geomlib_normalize_uvs(geomlib_uv_coords *uvs) {
    if (!uvs || !uvs->uvs || uvs->vertex_count == 0) return;

    /* Find bounding box */
    float min_u = FLT_MAX, min_v = FLT_MAX;
    float max_u = -FLT_MAX, max_v = -FLT_MAX;

    for (size_t i = 0; i < uvs->vertex_count; i++) {
        if (uvs->uvs[i].x < min_u) min_u = uvs->uvs[i].x;
        if (uvs->uvs[i].y < min_v) min_v = uvs->uvs[i].y;
        if (uvs->uvs[i].x > max_u) max_u = uvs->uvs[i].x;
        if (uvs->uvs[i].y > max_v) max_v = uvs->uvs[i].y;
    }

    float scale_u = (max_u > min_u) ? (1.0f / (max_u - min_u)) : 1.0f;
    float scale_v = (max_v > min_v) ? (1.0f / (max_v - min_v)) : 1.0f;

    /* Normalize to [0,1] */
    for (size_t i = 0; i < uvs->vertex_count; i++) {
        uvs->uvs[i].x = (uvs->uvs[i].x - min_u) * scale_u;
        uvs->uvs[i].y = (uvs->uvs[i].y - min_v) * scale_v;
    }
}

/* UV coordinates management */
geomlib_uv_coords geomlib_uv_coords_create(size_t vertex_count, size_t triangle_count) {
    geomlib_uv_coords uvs;
    uvs.uvs = calloc(vertex_count, sizeof(geomlib_vec2));
    uvs.vertex_count = vertex_count;
    uvs.distortion = triangle_count > 0 ? calloc(triangle_count, sizeof(float)) : NULL;
    uvs.triangle_count = triangle_count;
    return uvs;
}

void geomlib_uv_coords_free(geomlib_uv_coords *uvs) {
    if (!uvs) return;
    if (uvs->uvs) {
        free(uvs->uvs);
        uvs->uvs = NULL;
    }
    if (uvs->distortion) {
        free(uvs->distortion);
        uvs->distortion = NULL;
    }
    uvs->vertex_count = 0;
    uvs->triangle_count = 0;
}

/* Compute distortion metrics */
void geomlib_compute_param_distortion(
    const geomlib_trimesh *mesh,
    const geomlib_uv_coords *uvs
) {
    if (!mesh || !uvs || !uvs->distortion) return;

    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        /* 3D triangle */
        geomlib_vec3 p0 = mesh->vertices[i0];
        geomlib_vec3 p1 = mesh->vertices[i1];
        geomlib_vec3 p2 = mesh->vertices[i2];

        geomlib_vec3 e1_3d = geomlib_vec3_sub(p1, p0);
        geomlib_vec3 e2_3d = geomlib_vec3_sub(p2, p0);

        /* 2D triangle (UV) */
        geomlib_vec2 uv0 = uvs->uvs[i0];
        geomlib_vec2 uv1 = uvs->uvs[i1];
        geomlib_vec2 uv2 = uvs->uvs[i2];

        geomlib_vec2 e1_uv = geomlib_vec2_sub(uv1, uv0);
        geomlib_vec2 e2_uv = geomlib_vec2_sub(uv2, uv0);

        /* Area ratio */
        float area_3d = geomlib_vec3_length(geomlib_vec3_cross(e1_3d, e2_3d)) * 0.5f;
        float area_2d = fabsf(e1_uv.x * e2_uv.y - e1_uv.y * e2_uv.x) * 0.5f;

        /* Store distortion (area ratio) */
        uvs->distortion[t] = (area_3d > 1e-10f) ? (area_2d / area_3d) : 0.0f;
    }
}

/* Compute parameterization statistics */
void geomlib_compute_param_stats(
    const geomlib_trimesh *mesh,
    const geomlib_uv_coords *uvs,
    geomlib_param_stats *stats
) {
    if (!mesh || !uvs || !stats) return;

    memset(stats, 0, sizeof(geomlib_param_stats));

    float angle_error_sum = 0.0f;
    float area_error_sum = 0.0f;

    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        /* UV triangle */
        geomlib_vec2 uv0 = uvs->uvs[i0];
        geomlib_vec2 uv1 = uvs->uvs[i1];
        geomlib_vec2 uv2 = uvs->uvs[i2];

        geomlib_vec2 e1 = geomlib_vec2_sub(uv1, uv0);
        geomlib_vec2 e2 = geomlib_vec2_sub(uv2, uv0);

        /* Check for flips (negative area) */
        float signed_area = e1.x * e2.y - e1.y * e2.x;
        if (signed_area < 0.0f) {
            stats->has_flips = true;
            stats->flipped_count++;
        }

        /* Area distortion */
        if (uvs->distortion) {
            float ratio = uvs->distortion[t];
            area_error_sum += (ratio - 1.0f) * (ratio - 1.0f);
        }
    }

    stats->area_distortion = sqrtf(area_error_sum / mesh->triangle_count);
}
