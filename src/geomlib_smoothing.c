/*
 * geomlib - Mesh Smoothing Implementation
 */

#include "geomlib/geomlib_smoothing.h"
#include "geomlib/geomlib_healing.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/* Smooth mesh using specified method */
int geomlib_smooth_mesh(
    geomlib_trimesh *mesh,
    const geomlib_smoothing_params *params
) {
    if (!mesh || !params) return GEOMLIB_MESH_ERROR_INVALID;

    /* Build half-edge structure for neighbor queries */
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Identify boundary vertices */
    bool *is_boundary = NULL;
    if (params->preserve_boundaries) {
        is_boundary = calloc(mesh->vertex_count, sizeof(bool));

        for (size_t he = 0; he < hemesh.halfedge_count; he++) {
            if (hemesh.halfedges[he].twin == UINT32_MAX) {
                is_boundary[hemesh.halfedges[he].vertex] = true;
            }
        }
    }

    /* Allocate temporary vertex positions */
    geomlib_vec3 *new_positions = malloc(mesh->vertex_count * sizeof(geomlib_vec3));
    geomlib_vec3 *temp_positions = malloc(mesh->vertex_count * sizeof(geomlib_vec3));

    /* Initialize with current positions */
    memcpy(new_positions, mesh->vertices, mesh->vertex_count * sizeof(geomlib_vec3));

    /* Perform smoothing iterations */
    for (int iter = 0; iter < params->iterations; iter++) {
        /* Copy current positions to temp */
        memcpy(temp_positions, new_positions, mesh->vertex_count * sizeof(geomlib_vec3));

        /* Laplacian pass */
        for (size_t v = 0; v < mesh->vertex_count; v++) {
            /* Skip boundary vertices if preserving boundaries */
            if (params->preserve_boundaries && is_boundary && is_boundary[v]) {
                continue;
            }

            /* Compute average of neighbors */
            geomlib_vec3 sum = {0.0f, 0.0f, 0.0f};
            int neighbor_count = 0;

            /* Find one-ring neighbors */
            for (size_t he_idx = 0; he_idx < hemesh.halfedge_count; he_idx++) {
                if (hemesh.halfedges[he_idx].vertex != v) continue;

                uint32_t next = hemesh.halfedges[he_idx].next;
                if (next == UINT32_MAX) continue;

                uint32_t neighbor_idx = hemesh.halfedges[next].vertex;
                sum = geomlib_vec3_add(sum, temp_positions[neighbor_idx]);
                neighbor_count++;
            }

            if (neighbor_count > 0) {
                geomlib_vec3 avg = geomlib_vec3_div(sum, (float)neighbor_count);

                /* Update position based on method */
                geomlib_vec3 current = temp_positions[v];
                geomlib_vec3 delta = geomlib_vec3_sub(avg, current);

                if (params->method == GEOMLIB_SMOOTH_LAPLACIAN) {
                    /* Simple Laplacian: v' = v + lambda * (avg - v) */
                    new_positions[v] = geomlib_vec3_add(current,
                        geomlib_vec3_mul(delta, params->lambda));
                }
                else if (params->method == GEOMLIB_SMOOTH_TAUBIN) {
                    /* Taubin smoothing alternates between expansion and contraction */
                    float weight = (iter % 2 == 0) ? params->lambda : params->mu;
                    new_positions[v] = geomlib_vec3_add(current,
                        geomlib_vec3_mul(delta, weight));
                }
                else if (params->method == GEOMLIB_SMOOTH_HC) {
                    /* HC smoothing (simplified) */
                    geomlib_vec3 p = geomlib_vec3_add(current,
                        geomlib_vec3_mul(delta, params->alpha));
                    geomlib_vec3 b = geomlib_vec3_sub(p, current);
                    new_positions[v] = geomlib_vec3_sub(p,
                        geomlib_vec3_mul(b, params->beta));
                }
            }
        }
    }

    /* Update mesh vertices */
    memcpy(mesh->vertices, new_positions, mesh->vertex_count * sizeof(geomlib_vec3));

    /* Recompute normals */
    if (mesh->normals) {
        geomlib_trimesh_compute_normals(mesh);
    }

    /* Cleanup */
    free(new_positions);
    free(temp_positions);
    if (is_boundary) free(is_boundary);
    geomlib_halfedge_free(&hemesh);

    return GEOMLIB_MESH_SUCCESS;
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
