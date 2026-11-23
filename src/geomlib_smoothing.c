/*
 * geomlib - Mesh Smoothing Implementation
 */

#include "geomlib/geomlib_smoothing.h"
#include "geomlib/geomlib_feature.h"
#include "geomlib/geomlib_healing.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* Create default smoothing parameters */
geomlib_smoothing_params geomlib_smoothing_params_default(geomlib_smoothing_method method) {
    geomlib_smoothing_params params;
    params.method = method;
    params.lambda = 0.5f;
    params.iterations = 10;
    params.mu = -0.53f;
    params.alpha = 0.5f;
    params.beta = 0.5f;
    params.preserve_boundaries = true;
    params.preserve_features = false;
    params.feature_angle = 30.0f;
    return params;
}

static int compare_vertex_ids(const void *left, const void *right) {
    uint32_t a = *(const uint32_t*)left, b = *(const uint32_t*)right;
    return (a > b) - (a < b);
}

/* CSR one-ring: unique neighbors of every vertex. */
static int build_one_ring(const geomlib_trimesh *mesh, size_t **offsets_out,
                          uint32_t **neighbors_out, size_t **counts_out) {
    size_t vertices = mesh->vertex_count;
    size_t *counts = calloc(vertices, sizeof(*counts));
    size_t *offsets = calloc(vertices + 1, sizeof(*offsets));
    if (vertices > SIZE_MAX / (6 * sizeof(uint32_t)) ||
        mesh->triangle_count > SIZE_MAX / (6 * sizeof(uint32_t))) {
        free(counts);
        free(offsets);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    uint32_t *neighbors = malloc(mesh->triangle_count * 6 * sizeof(*neighbors));
    if (!counts || !offsets || !neighbors) {
        free(counts);
        free(offsets);
        free(neighbors);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        for (int j = 0; j < 3; j++) {
            uint32_t v = mesh->indices[t * 3 + j];
            if (v < vertices) counts[v] += 2;
        }
    }
    for (size_t i = 0; i < vertices; i++) offsets[i + 1] = offsets[i] + counts[i];
    size_t *cursor = malloc(vertices * sizeof(*cursor));
    if (!cursor) {
        free(counts);
        free(offsets);
        free(neighbors);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    memcpy(cursor, offsets, vertices * sizeof(*cursor));
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t v[3] = {mesh->indices[t * 3], mesh->indices[t * 3 + 1],
                         mesh->indices[t * 3 + 2]};
        for (int j = 0; j < 3; j++) {
            if (v[j] >= vertices) continue;
            neighbors[cursor[v[j]]++] = v[(j + 1) % 3];
            neighbors[cursor[v[j]]++] = v[(j + 2) % 3];
        }
    }
    free(cursor);
    for (size_t v = 0; v < vertices; v++) {
        size_t count = counts[v];
        uint32_t *start = neighbors + offsets[v];
        qsort(start, count, sizeof(*start), compare_vertex_ids);
        size_t unique = 0;
        for (size_t i = 0; i < count; i++)
            if (!unique || start[unique - 1] != start[i]) start[unique++] = start[i];
        counts[v] = unique;
    }
    *offsets_out = offsets;
    *neighbors_out = neighbors;
    *counts_out = counts;
    return GEOMLIB_MESH_SUCCESS;
}

static geomlib_vec3 average_position(const geomlib_vec3 *positions, const size_t *offsets,
                                     const uint32_t *neighbors, size_t count, uint32_t v) {
    size_t n = count;
    if (!n) return positions[v];
    geomlib_vec3d sum = {0};
    for (size_t i = 0; i < n; i++)
        sum = geomlib_vec3d_add(sum, (geomlib_vec3d){positions[neighbors[offsets[v] + i]].x,
                                                    positions[neighbors[offsets[v] + i]].y,
                                                    positions[neighbors[offsets[v] + i]].z});
    return (geomlib_vec3){(float)(sum.x / (double)n), (float)(sum.y / (double)n),
                          (float)(sum.z / (double)n)};
}

static void laplacian_pass(const geomlib_vec3 *in, geomlib_vec3 *out, float weight,
                           const bool *fixed, const size_t *offsets,
                           const uint32_t *neighbors, const size_t *counts,
                           size_t vertex_count) {
    for (uint32_t v = 0; v < vertex_count; v++) {
        if (fixed[v] || !counts[v]) {
            out[v] = in[v];
            continue;
        }
        geomlib_vec3 avg = average_position(in, offsets, neighbors, counts[v], v);
        out[v] = geomlib_vec3_add(in[v],
            geomlib_vec3_mul(geomlib_vec3_sub(avg, in[v]), weight));
    }
}

static bool *compute_fixed(const geomlib_trimesh *mesh, const geomlib_smoothing_params *params) {
    bool *fixed = calloc(mesh->vertex_count, sizeof(bool));
    if (!fixed) return NULL;
    if (params->preserve_boundaries) {
        geomlib_halfedge_mesh hemesh;
        if (geomlib_halfedge_build(&hemesh, mesh) == GEOMLIB_MESH_SUCCESS) {
            for (size_t he = 0; he < hemesh.halfedge_count; he++)
                if (hemesh.halfedges[he].twin == UINT32_MAX)
                    fixed[hemesh.halfedges[he].vertex] = true;
            geomlib_halfedge_free(&hemesh);
        }
    }
    if (params->preserve_features) {
        size_t count = 0;
        geomlib_feature_edge *edges =
            geomlib_detect_sharp_edges(mesh, params->feature_angle, &count);
        for (size_t i = 0; i < count; i++) {
            if (edges[i].v0 < mesh->vertex_count) fixed[edges[i].v0] = true;
            if (edges[i].v1 < mesh->vertex_count) fixed[edges[i].v1] = true;
        }
        free(edges);
    }
    return fixed;
}

/* Smooth mesh using specified method */
int geomlib_smooth_mesh(
    geomlib_trimesh *mesh,
    const geomlib_smoothing_params *params
) {
    if (!mesh || !params || params->iterations < 0 || !(params->lambda >= 0.0f) ||
        !(params->alpha >= 0.0f) || !(params->alpha <= 1.0f) || !isfinite(params->mu) ||
        !isfinite(params->lambda) || !isfinite(params->beta))
        return GEOMLIB_MESH_ERROR_INVALID;
    if (!params->iterations || !mesh->vertex_count || !mesh->triangle_count)
        return GEOMLIB_MESH_SUCCESS;

    size_t *offsets = NULL, *counts = NULL;
    uint32_t *neighbors = NULL;
    int result = build_one_ring(mesh, &offsets, &neighbors, &counts);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    bool *fixed = compute_fixed(mesh, params);
    geomlib_vec3 *new_positions = malloc(mesh->vertex_count * sizeof(*new_positions));
    geomlib_vec3 *temp = malloc(mesh->vertex_count * sizeof(*temp));
    geomlib_vec3 *p = NULL, *avg_p = NULL;
    if (!fixed || !new_positions || !temp) {
        result = GEOMLIB_MESH_ERROR_MEMORY;
        goto finish;
    }
    memcpy(new_positions, mesh->vertices, mesh->vertex_count * sizeof(*new_positions));

    if (params->method == GEOMLIB_SMOOTH_HC) {
        p = malloc(mesh->vertex_count * sizeof(*p));
        avg_p = malloc(mesh->vertex_count * sizeof(*avg_p));
        if (!p || !avg_p) {
            result = GEOMLIB_MESH_ERROR_MEMORY;
            goto finish;
        }
    }

    for (int iter = 0; iter < params->iterations; iter++) {
        if (params->method == GEOMLIB_SMOOTH_TAUBIN) {
            memcpy(temp, new_positions, mesh->vertex_count * sizeof(*temp));
            laplacian_pass(temp, new_positions, params->lambda, fixed, offsets, neighbors,
                           counts, mesh->vertex_count);
            memcpy(temp, new_positions, mesh->vertex_count * sizeof(*temp));
            laplacian_pass(temp, new_positions, params->mu, fixed, offsets, neighbors,
                           counts, mesh->vertex_count);
        } else if (params->method == GEOMLIB_SMOOTH_HC) {
            memcpy(temp, new_positions, mesh->vertex_count * sizeof(*temp));
            for (uint32_t v = 0; v < mesh->vertex_count; v++) {
                if (fixed[v] || !counts[v]) {
                    p[v] = temp[v];
                    continue;
                }
                geomlib_vec3 avg = average_position(temp, offsets, neighbors, counts[v], v);
                p[v] = geomlib_vec3_add(temp[v],
                    geomlib_vec3_mul(geomlib_vec3_sub(avg, temp[v]), params->alpha));
            }
            for (uint32_t v = 0; v < mesh->vertex_count; v++)
                avg_p[v] = counts[v] ? average_position(p, offsets, neighbors, counts[v], v) : p[v];
            for (uint32_t v = 0; v < mesh->vertex_count; v++) {
                if (fixed[v] || !counts[v]) {
                    new_positions[v] = p[v];
                    continue;
                }
                geomlib_vec3 b = geomlib_vec3_sub(p[v], geomlib_vec3_add(
                    geomlib_vec3_mul(p[v], params->alpha),
                    geomlib_vec3_mul(avg_p[v], 1.0f - params->alpha)));
                new_positions[v] = geomlib_vec3_sub(p[v], geomlib_vec3_mul(b, params->beta));
            }
        } else {
            memcpy(temp, new_positions, mesh->vertex_count * sizeof(*temp));
            laplacian_pass(temp, new_positions, params->lambda, fixed, offsets, neighbors,
                           counts, mesh->vertex_count);
        }
    }

    memcpy(mesh->vertices, new_positions, mesh->vertex_count * sizeof(*new_positions));
    if (mesh->normals) geomlib_trimesh_compute_normals(mesh);

finish:
    free(offsets);
    free(neighbors);
    free(counts);
    free(fixed);
    free(new_positions);
    free(temp);
    free(p);
    free(avg_p);
    return result;
}

/* Simplified interfaces */

int geomlib_smooth_laplacian(
    geomlib_trimesh *mesh,
    float lambda,
    int iterations
) {
    geomlib_smoothing_params params = geomlib_smoothing_params_default(GEOMLIB_SMOOTH_LAPLACIAN);
    params.lambda = lambda;
    params.iterations = iterations;
    return geomlib_smooth_mesh(mesh, &params);
}

int geomlib_smooth_taubin(
    geomlib_trimesh *mesh,
    float lambda,
    float mu,
    int iterations
) {
    geomlib_smoothing_params params = geomlib_smoothing_params_default(GEOMLIB_SMOOTH_TAUBIN);
    params.lambda = lambda;
    params.mu = mu;
    params.iterations = iterations;
    return geomlib_smooth_mesh(mesh, &params);
}

int geomlib_smooth_normals(
    geomlib_trimesh *mesh,
    int iterations
) {
    if (!mesh || !mesh->normals) return GEOMLIB_MESH_ERROR_INVALID;

    /* Build half-edge structure */
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    geomlib_vec3 *new_normals = malloc(mesh->vertex_count * sizeof(geomlib_vec3));
    memcpy(new_normals, mesh->normals, mesh->vertex_count * sizeof(geomlib_vec3));

    /* Smooth normals */
    for (int iter = 0; iter < iterations; iter++) {
        geomlib_vec3 *temp_normals = malloc(mesh->vertex_count * sizeof(geomlib_vec3));
        memcpy(temp_normals, new_normals, mesh->vertex_count * sizeof(geomlib_vec3));

        for (size_t v = 0; v < mesh->vertex_count; v++) {
            geomlib_vec3 sum = {0.0f, 0.0f, 0.0f};
            int count = 0;

            /* Average with neighbors */
            for (size_t he_idx = 0; he_idx < hemesh.halfedge_count; he_idx++) {
                if (hemesh.halfedges[he_idx].vertex != v) continue;

                uint32_t next = hemesh.halfedges[he_idx].next;
                if (next == UINT32_MAX) continue;

                uint32_t neighbor_idx = hemesh.halfedges[next].vertex;
                sum = geomlib_vec3_add(sum, temp_normals[neighbor_idx]);
                count++;
            }

            if (count > 0) {
                new_normals[v] = geomlib_vec3_normalize(
                    geomlib_vec3_div(sum, (float)count));
            }
        }

        free(temp_normals);
    }

    memcpy(mesh->normals, new_normals, mesh->vertex_count * sizeof(geomlib_vec3));
    free(new_normals);
    geomlib_halfedge_free(&hemesh);

    return GEOMLIB_MESH_SUCCESS;
}
