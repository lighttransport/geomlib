/*
 * geomlib - QEM Decimation Implementation
 * Simplified version of Garland-Heckbert algorithm
 */

#include "geomlib/geomlib_decimation.h"
#include "geomlib/geomlib_healing.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* Create default decimation parameters */
geomlib_decimation_params geomlib_decimation_params_default(float target_reduction) {
    geomlib_decimation_params params;
    params.target_reduction = target_reduction;
    params.target_vertex_count = 0;
    params.target_triangle_count = 0;
    params.max_error = 0.0f;  /* No limit */
    params.preserve_boundaries = true;
    params.preserve_topology = true;
    params.boundary_weight = 1000.0f;
    return params;
}

/* Initialize a quadric from a plane */
geomlib_quadric geomlib_quadric_from_plane(float a, float b, float c, float d) {
    geomlib_quadric q;
    q.a2 = a * a;
    q.ab = a * b;
    q.ac = a * c;
    q.ad = a * d;
    q.b2 = b * b;
    q.bc = b * c;
    q.bd = b * d;
    q.c2 = c * c;
    q.cd = c * d;
    q.d2 = d * d;
    return q;
}

/* Add two quadrics */
geomlib_quadric geomlib_quadric_add(geomlib_quadric q1, geomlib_quadric q2) {
    geomlib_quadric result;
    result.a2 = q1.a2 + q2.a2;
    result.ab = q1.ab + q2.ab;
    result.ac = q1.ac + q2.ac;
    result.ad = q1.ad + q2.ad;
    result.b2 = q1.b2 + q2.b2;
    result.bc = q1.bc + q2.bc;
    result.bd = q1.bd + q2.bd;
    result.c2 = q1.c2 + q2.c2;
    result.cd = q1.cd + q2.cd;
    result.d2 = q1.d2 + q2.d2;
    return result;
}

/* Compute quadric error for a vertex position */
double geomlib_quadric_error(geomlib_quadric q, geomlib_vec3 v) {
    double x = v.x, y = v.y, z = v.z;

    double error =
        q.a2 * x * x + 2.0 * q.ab * x * y + 2.0 * q.ac * x * z + 2.0 * q.ad * x +
        q.b2 * y * y + 2.0 * q.bc * y * z + 2.0 * q.bd * y +
        q.c2 * z * z + 2.0 * q.cd * z +
        q.d2;

    return error;
}

/* Simplified QEM decimation using greedy edge collapse */
int geomlib_decimate_qem(
    geomlib_trimesh *mesh,
    const geomlib_decimation_params *params,
    geomlib_decimation_stats *stats
) {
    if (!mesh || !params) return GEOMLIB_MESH_ERROR_INVALID;

    /* Determine target vertex count */
    size_t target_verts = mesh->vertex_count;

    if (params->target_vertex_count > 0) {
        target_verts = params->target_vertex_count;
    } else if (params->target_triangle_count > 0) {
        /* Euler formula: V - E + F = 2 for closed mesh */
        /* For triangle mesh: E ≈ 1.5 * F, so V ≈ F/2 + 2 */
        target_verts = params->target_triangle_count / 2 + 2;
    } else {
        target_verts = (size_t)((1.0f - params->target_reduction) * mesh->vertex_count);
    }

    if (target_verts >= mesh->vertex_count) {
        return GEOMLIB_MESH_SUCCESS;  /* Nothing to do */
    }

    /* Initialize statistics */
    geomlib_decimation_stats local_stats = {0};
    local_stats.original_vertices = mesh->vertex_count;
    local_stats.original_triangles = mesh->triangle_count;

    /* Build half-edge structure */
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Compute quadric for each vertex */
    geomlib_quadric *vertex_quadrics = calloc(mesh->vertex_count, sizeof(geomlib_quadric));
    if (!vertex_quadrics) {
        geomlib_halfedge_free(&hemesh);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    /* Initialize quadrics from incident triangles */
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        geomlib_vec3 p0 = mesh->vertices[i0];
        geomlib_vec3 p1 = mesh->vertices[i1];
        geomlib_vec3 p2 = mesh->vertices[i2];

        /* Compute triangle normal */
        geomlib_vec3 e1 = geomlib_vec3_sub(p1, p0);
        geomlib_vec3 e2 = geomlib_vec3_sub(p2, p0);
        geomlib_vec3 normal = geomlib_vec3_normalize(geomlib_vec3_cross(e1, e2));

        /* Compute plane equation: ax + by + cz + d = 0 */
        float a = normal.x;
        float b = normal.y;
        float c = normal.z;
        float d = -(a * p0.x + b * p0.y + c * p0.z);

        /* Create quadric for this plane */
        geomlib_quadric tri_quad = geomlib_quadric_from_plane(a, b, c, d);

        /* Add to vertex quadrics */
        vertex_quadrics[i0] = geomlib_quadric_add(vertex_quadrics[i0], tri_quad);
        vertex_quadrics[i1] = geomlib_quadric_add(vertex_quadrics[i1], tri_quad);
        vertex_quadrics[i2] = geomlib_quadric_add(vertex_quadrics[i2], tri_quad);
    }

    /* Mark vertices to keep (simple greedy approach) */
    /* In a full implementation, we'd use a priority queue of edge collapses */
    bool *vertex_removed = calloc(mesh->vertex_count, sizeof(bool));
    size_t removed_count = 0;
    size_t vertices_to_remove = mesh->vertex_count - target_verts;

    /* Greedy removal: remove vertices with lowest quadric error */
    while (removed_count < vertices_to_remove) {
        /* Find vertex with minimum error when moved to its position */
        uint32_t best_vertex = UINT32_MAX;
        double min_error = DBL_MAX;

        for (size_t v = 0; v < mesh->vertex_count; v++) {
            if (vertex_removed[v]) continue;

            /* Skip boundary vertices if preserving boundaries */
            if (params->preserve_boundaries) {
                /* Check if this is a boundary vertex */
                bool is_boundary = false;
                for (size_t he = 0; he < hemesh.halfedge_count; he++) {
                    if (hemesh.halfedges[he].vertex == v && hemesh.halfedges[he].twin == UINT32_MAX) {
                        is_boundary = true;
                        break;
                    }
                }
                if (is_boundary) continue;
            }

            double error = geomlib_quadric_error(vertex_quadrics[v], mesh->vertices[v]);

            if (error < min_error) {
                min_error = error;
                best_vertex = v;
            }
        }

        if (best_vertex == UINT32_MAX) break;  /* No more removable vertices */

        /* Check error threshold */
        if (params->max_error > 0.0f && min_error > params->max_error) {
            break;  /* Error threshold exceeded */
        }

        /* Mark vertex for removal */
        vertex_removed[best_vertex] = true;
        removed_count++;

        /* Update statistics */
        if (min_error > local_stats.max_error) {
            local_stats.max_error = (float)min_error;
        }
        local_stats.mean_error += (float)min_error;
    }

    if (removed_count > 0) {
        local_stats.mean_error /= removed_count;
    }

    /* Build new mesh without removed vertices */
    uint32_t *vertex_map = malloc(mesh->vertex_count * sizeof(uint32_t));
    uint32_t new_vertex_count = 0;

    for (size_t v = 0; v < mesh->vertex_count; v++) {
        if (!vertex_removed[v]) {
            vertex_map[v] = new_vertex_count++;
        } else {
            vertex_map[v] = UINT32_MAX;
        }
    }

    /* Create new vertex array */
    geomlib_vec3 *new_vertices = malloc(new_vertex_count * sizeof(geomlib_vec3));
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        if (!vertex_removed[v]) {
            new_vertices[vertex_map[v]] = mesh->vertices[v];
        }
    }

    /* Create new triangle array (skip triangles with removed vertices) */
    uint32_t *new_indices = malloc(mesh->triangle_count * 3 * sizeof(uint32_t));
    size_t new_triangle_count = 0;

    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3 + 0];
        uint32_t i1 = mesh->indices[t * 3 + 1];
        uint32_t i2 = mesh->indices[t * 3 + 2];

        /* Skip if any vertex was removed */
        if (vertex_removed[i0] || vertex_removed[i1] || vertex_removed[i2]) {
            continue;
        }

        /* Remap indices */
        new_indices[new_triangle_count * 3 + 0] = vertex_map[i0];
        new_indices[new_triangle_count * 3 + 1] = vertex_map[i1];
        new_indices[new_triangle_count * 3 + 2] = vertex_map[i2];
        new_triangle_count++;
    }

    /* Update mesh */
    free(mesh->vertices);
    free(mesh->indices);
    if (mesh->normals) {
        free(mesh->normals);
        mesh->normals = NULL;
    }

    mesh->vertices = new_vertices;
    mesh->indices = new_indices;
    mesh->vertex_count = new_vertex_count;
    mesh->triangle_count = new_triangle_count;

    /* Update capacity */
    mesh->vertex_capacity = new_vertex_count;
    mesh->triangle_capacity = new_triangle_count;

    /* Cleanup */
    free(vertex_quadrics);
    free(vertex_removed);
    free(vertex_map);
    geomlib_halfedge_free(&hemesh);

    /* Finalize statistics */
    local_stats.final_vertices = new_vertex_count;
    local_stats.final_triangles = new_triangle_count;
    local_stats.edges_collapsed = removed_count;
    local_stats.reduction_ratio = 1.0f - ((float)new_vertex_count / local_stats.original_vertices);

    if (stats) {
        *stats = local_stats;
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* Simplified interfaces */
int geomlib_decimate_to_ratio(geomlib_trimesh *mesh, float target_ratio) {
    geomlib_decimation_params params = geomlib_decimation_params_default(1.0f - target_ratio);
    return geomlib_decimate_qem(mesh, &params, NULL);
}

int geomlib_decimate_to_vertex_count(geomlib_trimesh *mesh, size_t target_vertices) {
    geomlib_decimation_params params = geomlib_decimation_params_default(0.0f);
    params.target_vertex_count = target_vertices;
    return geomlib_decimate_qem(mesh, &params, NULL);
}

int geomlib_decimate_to_triangle_count(geomlib_trimesh *mesh, size_t target_triangles) {
    geomlib_decimation_params params = geomlib_decimation_params_default(0.0f);
    params.target_triangle_count = target_triangles;
    return geomlib_decimate_qem(mesh, &params, NULL);
}
