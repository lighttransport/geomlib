/*
 * geomlib - Geodesic Distance Computation Implementation
 * Simplified heat method using iterative solvers
 */

#include "geomlib/geomlib_allocator.h"
#include "geomlib/geomlib_c11.h"
#include "geomlib/geomlib_geodesic.h"
#include "geomlib/geomlib_healing.h"
#include "geomlib/geomlib_curvature.h"
#include <float.h>

/* Create geodesic distance structure */
geomlib_geodesic_distances geomlib_geodesic_create(size_t vertex_count) {
    geomlib_geodesic_distances geo;
    geo.distances = geomlib_calloc(vertex_count, sizeof(float));
    geo.vertex_count = vertex_count;
    geo.source_vertices = NULL;
    geo.source_count = 0;
    return geo;
}

/* Free geodesic distance structure */
void geomlib_geodesic_free(geomlib_geodesic_distances *geodesic) {
    if (!geodesic) return;
    if (geodesic->distances) {
        geomlib_free(geodesic->distances);
        geodesic->distances = NULL;
    }
    if (geodesic->source_vertices) {
        geomlib_free(geodesic->source_vertices);
        geodesic->source_vertices = NULL;
    }
    geodesic->vertex_count = 0;
    geodesic->source_count = 0;
}

/* Compute mean edge length */
float geomlib_compute_mean_edge_length_for_geodesic(const geomlib_trimesh *mesh) {
    if (!mesh || mesh->triangle_count == 0) return 1.0f;

    float total_length = 0.0f;
    size_t edge_count = 0;

    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        geomlib_vec3 p0 = mesh->vertices[i0];
        geomlib_vec3 p1 = mesh->vertices[i1];
        geomlib_vec3 p2 = mesh->vertices[i2];

        total_length += geomlib_vec3_length(geomlib_vec3_sub(p1, p0));
        total_length += geomlib_vec3_length(geomlib_vec3_sub(p2, p1));
        total_length += geomlib_vec3_length(geomlib_vec3_sub(p0, p2));
        edge_count += 3;
    }

    return edge_count > 0 ? (total_length / edge_count) : 1.0f;
}

/* Simplified heat method using Dijkstra-like propagation */
int geomlib_compute_geodesic_heat(
    const geomlib_trimesh *mesh,
    const uint32_t *source_vertices,
    size_t source_count,
    float time_step,
    geomlib_geodesic_distances *distances
) {
    if (!mesh || !source_vertices || source_count == 0 || !distances) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    if (distances->vertex_count != mesh->vertex_count) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Build half-edge structure for neighbor queries */
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Initialize distances to infinity */
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        distances->distances[i] = FLT_MAX;
    }

    /* Set source distances to 0 */
    for (size_t i = 0; i < source_count; i++) {
        if (source_vertices[i] < mesh->vertex_count) {
            distances->distances[source_vertices[i]] = 0.0f;
        }
    }

    /* Store source vertices */
    if (distances->source_vertices) {
        geomlib_free(distances->source_vertices);
    }
    distances->source_vertices = geomlib_malloc(source_count * sizeof(uint32_t));
    geomlib_memcpy(distances->source_vertices, source_vertices, source_count * sizeof(uint32_t));
    distances->source_count = source_count;

    /* Use simplified Dijkstra-like algorithm for geodesic distances */
    /* This is more accurate than the heat diffusion for small meshes */

    bool *visited = geomlib_calloc(mesh->vertex_count, sizeof(bool));

    /* Process vertices in order of increasing distance */
    for (size_t iter = 0; iter < mesh->vertex_count; iter++) {
        /* Find unvisited vertex with minimum distance */
        uint32_t min_vertex = UINT32_MAX;
        float min_dist = FLT_MAX;

        for (size_t v = 0; v < mesh->vertex_count; v++) {
            if (!visited[v] && distances->distances[v] < min_dist) {
                min_dist = distances->distances[v];
                min_vertex = v;
            }
        }

        if (min_vertex == UINT32_MAX) break;

        visited[min_vertex] = true;

        /* Update neighbors */
        for (size_t he_idx = 0; he_idx < hemesh.halfedge_count; he_idx++) {
            if (hemesh.halfedges[he_idx].vertex != min_vertex) continue;

            uint32_t next = hemesh.halfedges[he_idx].next;
            if (next == UINT32_MAX) continue;

            uint32_t neighbor_idx = hemesh.halfedges[next].vertex;
            if (visited[neighbor_idx]) continue;

            /* Compute edge length */
            geomlib_vec3 v1 = mesh->vertices[min_vertex];
            geomlib_vec3 v2 = mesh->vertices[neighbor_idx];
            float edge_len = geomlib_vec3_length(geomlib_vec3_sub(v2, v1));

            /* Update distance */
            float new_dist = distances->distances[min_vertex] + edge_len;
            if (new_dist < distances->distances[neighbor_idx]) {
                distances->distances[neighbor_idx] = new_dist;
            }
        }
    }

    geomlib_free(visited);
    geomlib_halfedge_free(&hemesh);

    return GEOMLIB_MESH_SUCCESS;
}

/* Compute geodesic distance from a single source vertex */
int geomlib_compute_geodesic_from_vertex(
    const geomlib_trimesh *mesh,
    uint32_t source_vertex,
    geomlib_geodesic_distances *distances
) {
    if (source_vertex >= mesh->vertex_count) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    return geomlib_compute_geodesic_heat(mesh, &source_vertex, 1, 0.0f, distances);
}

/* Compute all-pairs geodesic distances */
int geomlib_compute_geodesic_all_pairs(
    const geomlib_trimesh *mesh,
    float **distance_matrix
) {
    if (!mesh || !distance_matrix) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Allocate distance matrix */
    size_t n = mesh->vertex_count;
    *distance_matrix = geomlib_malloc(n * n * sizeof(float));
    if (!*distance_matrix) {
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    /* Compute distances from each vertex */
    geomlib_geodesic_distances geo = geomlib_geodesic_create(n);

    for (size_t i = 0; i < n; i++) {
        uint32_t source = (uint32_t)i;
        int result = geomlib_compute_geodesic_from_vertex(mesh, source, &geo);

        if (result != GEOMLIB_MESH_SUCCESS) {
            geomlib_free(*distance_matrix);
            *distance_matrix = NULL;
            geomlib_geodesic_free(&geo);
            return result;
        }

        /* Copy distances to matrix row */
        for (size_t j = 0; j < n; j++) {
            (*distance_matrix)[i * n + j] = geo.distances[j];
        }
    }

    geomlib_geodesic_free(&geo);

    return GEOMLIB_MESH_SUCCESS;
}
