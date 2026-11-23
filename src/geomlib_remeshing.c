/*
 * geomlib - Remeshing Implementation
 * Based on Botsch-Kobbelt incremental remeshing (SGP 2004)
 */

#include "geomlib/geomlib_remeshing.h"
#include "geomlib/geomlib_healing.h"
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Helper: Compute edge length */
static float edge_length(const geomlib_trimesh *mesh, uint32_t v0, uint32_t v1) {
    geomlib_vec3 p0 = mesh->vertices[v0];
    geomlib_vec3 p1 = mesh->vertices[v1];
    geomlib_vec3 diff = geomlib_vec3_sub(p1, p0);
    return geomlib_vec3_length(diff);
}

/* Helper: Compute triangle normal */
static geomlib_vec3 triangle_normal(const geomlib_trimesh *mesh, uint32_t tri_idx) {
    uint32_t i0 = mesh->indices[tri_idx * 3 + 0];
    uint32_t i1 = mesh->indices[tri_idx * 3 + 1];
    uint32_t i2 = mesh->indices[tri_idx * 3 + 2];

    geomlib_vec3 p0 = mesh->vertices[i0];
    geomlib_vec3 p1 = mesh->vertices[i1];
    geomlib_vec3 p2 = mesh->vertices[i2];

    geomlib_vec3 e1 = geomlib_vec3_sub(p1, p0);
    geomlib_vec3 e2 = geomlib_vec3_sub(p2, p0);
    geomlib_vec3 n = geomlib_vec3_cross(e1, e2);
    return geomlib_vec3_normalize(n);
}

/* Helper: Compute vertex normal from face normals */
static geomlib_vec3 compute_vertex_normal(
    const geomlib_trimesh *mesh,
    const geomlib_halfedge_mesh *hemesh,
    uint32_t vertex_idx
) {
    geomlib_vec3 normal = geomlib_vec3_create(0, 0, 0);
    int face_count = 0;

    /* Find a half-edge starting from this vertex */
    uint32_t he_start = UINT32_MAX;
    for (size_t i = 0; i < hemesh->halfedge_count; i++) {
        if (hemesh->halfedges[i].vertex == vertex_idx) {
            he_start = i;
            break;
        }
    }

    if (he_start == UINT32_MAX) return geomlib_vec3_create(0, 1, 0);

    /* Traverse one-ring */
    uint32_t he_current = he_start;
    do {
        uint32_t face = hemesh->halfedges[he_current].face;
        if (face != UINT32_MAX) {
            geomlib_vec3 fn = triangle_normal(mesh, face);
            normal = geomlib_vec3_add(normal, fn);
            face_count++;
        }

        /* Move to twin->next */
        uint32_t twin = hemesh->halfedges[he_current].twin;
        if (twin == UINT32_MAX) break;
        he_current = hemesh->halfedges[twin].next;

        if (he_current == he_start) break;
    } while (he_current != UINT32_MAX && face_count < 100);

    if (face_count > 0) {
        normal = geomlib_vec3_mul(normal, 1.0f / face_count);
    }

    return geomlib_vec3_normalize(normal);
}

/* Helper: Check if edge collapse would create degenerate triangles */
static bool is_collapse_valid(
    const geomlib_trimesh *mesh,
    const geomlib_halfedge_mesh *hemesh,
    uint32_t v0,
    uint32_t v1,
    geomlib_vec3 new_pos
) {
    /* Check all faces incident to v0 */
    for (size_t i = 0; i < hemesh->halfedge_count; i++) {
        if (hemesh->halfedges[i].vertex == v0) {
            uint32_t face = hemesh->halfedges[i].face;
            if (face == UINT32_MAX) continue;

            /* Get triangle vertices */
            uint32_t i0 = mesh->indices[face * 3 + 0];
            uint32_t i1 = mesh->indices[face * 3 + 1];
            uint32_t i2 = mesh->indices[face * 3 + 2];

            /* Skip faces shared with v1 (they will be removed) */
            if (i0 == v1 || i1 == v1 || i2 == v1) continue;

            /* Compute new triangle with v0 replaced by new_pos */
            geomlib_vec3 p0 = (i0 == v0) ? new_pos : mesh->vertices[i0];
            geomlib_vec3 p1 = (i1 == v0) ? new_pos : mesh->vertices[i1];
            geomlib_vec3 p2 = (i2 == v0) ? new_pos : mesh->vertices[i2];

            /* Check for degeneracy */
            geomlib_vec3 e1 = geomlib_vec3_sub(p1, p0);
            geomlib_vec3 e2 = geomlib_vec3_sub(p2, p0);
            geomlib_vec3 n = geomlib_vec3_cross(e1, e2);
            float area = geomlib_vec3_length(n);

            if (area < 1e-10f) return false;
        }
    }

    return true;
}

/* Helper: Get vertex valence (number of incident edges) */
static int get_vertex_valence(const geomlib_halfedge_mesh *hemesh, uint32_t vertex_idx) {
    int valence = 0;

    /* Find a half-edge starting from this vertex */
    uint32_t he_start = UINT32_MAX;
    for (size_t i = 0; i < hemesh->halfedge_count; i++) {
        if (hemesh->halfedges[i].vertex == vertex_idx) {
            he_start = i;
            break;
        }
    }

    if (he_start == UINT32_MAX) return 0;

    /* Count edges in one-ring */
    uint32_t he_current = he_start;
    do {
        valence++;

        /* Move to twin->next */
        uint32_t twin = hemesh->halfedges[he_current].twin;
        if (twin == UINT32_MAX) {
            /* Boundary vertex - count incoming edge too */
            return valence + 1;
        }
        he_current = hemesh->halfedges[twin].next;

        if (he_current == he_start) break;
    } while (he_current != UINT32_MAX && valence < 100);

    return valence;
}


/* Edge split operation */
static bool split_edge(
    geomlib_trimesh *mesh,
    geomlib_halfedge_mesh *hemesh,
    uint32_t v0,
    uint32_t v1
) {
    /* Find the half-edge */
    uint32_t he_idx = UINT32_MAX;
    for (size_t i = 0; i < hemesh->halfedge_count; i++) {
        uint32_t next = hemesh->halfedges[i].next;
        if (next != UINT32_MAX &&
            hemesh->halfedges[i].vertex == v0 &&
            hemesh->halfedges[next].vertex == v1) {
            he_idx = i;
            break;
        }
    }

    if (he_idx == UINT32_MAX) return false;

    /* Compute midpoint */
    geomlib_vec3 p0 = mesh->vertices[v0];
    geomlib_vec3 p1 = mesh->vertices[v1];
    geomlib_vec3 mid = geomlib_vec3_mul(geomlib_vec3_add(p0, p1), 0.5f);

    /* Add new vertex */
    uint32_t v_new = mesh->vertex_count;
    geomlib_trimesh_add_vertex(mesh, mid, NULL, NULL);

    /* Get affected faces */
    uint32_t face0 = hemesh->halfedges[he_idx].face;
    uint32_t twin_idx = hemesh->halfedges[he_idx].twin;
    uint32_t face1 = (twin_idx != UINT32_MAX) ? hemesh->halfedges[twin_idx].face : UINT32_MAX;

    /* Split face0 */
    if (face0 != UINT32_MAX) {
        uint32_t i0 = mesh->indices[face0 * 3 + 0];
        uint32_t i1 = mesh->indices[face0 * 3 + 1];
        uint32_t i2 = mesh->indices[face0 * 3 + 2];

        /* Find the opposite vertex */

        /* Modify face0: v0, v_new, v_opp */
        /* Add new face: v_new, v1, v_opp */
        if (i0 == v0 && i1 == v1) {
            mesh->indices[face0 * 3 + 1] = v_new;
            geomlib_trimesh_add_triangle(mesh, v_new, v1, i2);
        } else if (i1 == v0 && i2 == v1) {
            mesh->indices[face0 * 3 + 2] = v_new;
            geomlib_trimesh_add_triangle(mesh, v_new, v1, i0);
        } else if (i2 == v0 && i0 == v1) {
            mesh->indices[face0 * 3 + 0] = v_new;
            geomlib_trimesh_add_triangle(mesh, v_new, v1, i1);
        }
    }

    /* Split face1 if it exists */
    if (face1 != UINT32_MAX) {
        uint32_t i0 = mesh->indices[face1 * 3 + 0];
        uint32_t i1 = mesh->indices[face1 * 3 + 1];
        uint32_t i2 = mesh->indices[face1 * 3 + 2];


        if (i0 == v1 && i1 == v0) {
            mesh->indices[face1 * 3 + 1] = v_new;
            geomlib_trimesh_add_triangle(mesh, v_new, v0, i2);
        } else if (i1 == v1 && i2 == v0) {
            mesh->indices[face1 * 3 + 2] = v_new;
            geomlib_trimesh_add_triangle(mesh, v_new, v0, i0);
        } else if (i2 == v1 && i0 == v0) {
            mesh->indices[face1 * 3 + 0] = v_new;
            geomlib_trimesh_add_triangle(mesh, v_new, v0, i1);
        }
    }

    /* Rebuild half-edge structure */
    geomlib_halfedge_build(hemesh, mesh);

    return true;
}

/* Tangential smoothing */
static void tangential_smoothing(
    geomlib_trimesh *mesh,
    const geomlib_halfedge_mesh *hemesh,
    bool *is_feature,
    bool preserve_boundaries
) {
    geomlib_vec3 *new_positions = malloc(mesh->vertex_count * sizeof(geomlib_vec3));
    memcpy(new_positions, mesh->vertices, mesh->vertex_count * sizeof(geomlib_vec3));

    for (size_t v = 0; v < mesh->vertex_count; v++) {
        /* Skip feature vertices */
        if (is_feature && is_feature[v]) continue;

        /* Find one-ring neighbors */
        uint32_t he_start = UINT32_MAX;
        for (size_t i = 0; i < hemesh->halfedge_count; i++) {
            if (hemesh->halfedges[i].vertex == v) {
                he_start = i;
                break;
            }
        }

        if (he_start == UINT32_MAX) continue;

        /* Check if boundary vertex */
        bool is_boundary = false;
        uint32_t he_current = he_start;
        do {
            uint32_t twin = hemesh->halfedges[he_current].twin;
            if (twin == UINT32_MAX) {
                is_boundary = true;
                break;
            }
            he_current = hemesh->halfedges[twin].next;
        } while (he_current != he_start && he_current != UINT32_MAX);

        if (is_boundary && preserve_boundaries) continue;

        /* Compute centroid of neighbors */
        geomlib_vec3 centroid = geomlib_vec3_create(0, 0, 0);
        int neighbor_count = 0;

        he_current = he_start;
        do {
            uint32_t next = hemesh->halfedges[he_current].next;
            if (next != UINT32_MAX) {
                uint32_t neighbor = hemesh->halfedges[next].vertex;
                centroid = geomlib_vec3_add(centroid, mesh->vertices[neighbor]);
                neighbor_count++;
            }

            uint32_t twin = hemesh->halfedges[he_current].twin;
            if (twin == UINT32_MAX) break;
            he_current = hemesh->halfedges[twin].next;

            if (he_current == he_start) break;
        } while (he_current != UINT32_MAX && neighbor_count < 100);

        if (neighbor_count == 0) continue;

        centroid = geomlib_vec3_mul(centroid, 1.0f / neighbor_count);

        /* Compute vertex normal */
        geomlib_vec3 normal = compute_vertex_normal(mesh, hemesh, v);

        /* Project displacement to tangent plane */
        geomlib_vec3 displacement = geomlib_vec3_sub(centroid, mesh->vertices[v]);
        float proj = geomlib_vec3_dot(displacement, normal);
        geomlib_vec3 tangent_disp = geomlib_vec3_sub(displacement,
                                                      geomlib_vec3_mul(normal, proj));

        /* Apply smoothing */
        new_positions[v] = geomlib_vec3_add(mesh->vertices[v], tangent_disp);
    }

    /* Update positions */
    memcpy(mesh->vertices, new_positions, mesh->vertex_count * sizeof(geomlib_vec3));
    free(new_positions);
}

/* Detect feature edges based on dihedral angle */
static bool* detect_features(
    const geomlib_trimesh *mesh,
    const geomlib_halfedge_mesh *hemesh,
    float feature_angle_degrees
) {
    bool *is_feature = calloc(mesh->vertex_count, sizeof(bool));
    float cos_threshold = cosf(feature_angle_degrees * M_PI / 180.0f);

    for (size_t he_idx = 0; he_idx < hemesh->halfedge_count; he_idx++) {
        uint32_t twin = hemesh->halfedges[he_idx].twin;
        if (twin == UINT32_MAX || twin < he_idx) continue; /* Skip boundaries and duplicates */

        uint32_t face0 = hemesh->halfedges[he_idx].face;
        uint32_t face1 = hemesh->halfedges[twin].face;

        if (face0 == UINT32_MAX || face1 == UINT32_MAX) continue;

        /* Compute face normals */
        geomlib_vec3 n0 = triangle_normal(mesh, face0);
        geomlib_vec3 n1 = triangle_normal(mesh, face1);

        float cos_angle = geomlib_vec3_dot(n0, n1);

        if (cos_angle < cos_threshold) {
            /* Mark vertices of this edge as feature */
            uint32_t v0 = hemesh->halfedges[he_idx].vertex;
            uint32_t next = hemesh->halfedges[he_idx].next;
            uint32_t v1 = (next != UINT32_MAX) ? hemesh->halfedges[next].vertex : UINT32_MAX;

            if (v0 != UINT32_MAX) is_feature[v0] = true;
            if (v1 != UINT32_MAX) is_feature[v1] = true;
        }
    }

    return is_feature;
}

/* Compute curvature-based target edge length */
static float compute_adaptive_edge_length(
    const geomlib_trimesh *mesh,
    const geomlib_halfedge_mesh *hemesh,
    uint32_t vertex_idx,
    float min_length,
    float max_length
) {
    /* Estimate curvature from vertex normal variation */
    geomlib_vec3 normal = compute_vertex_normal(mesh, hemesh, vertex_idx);

    /* Find one-ring */
    uint32_t he_start = UINT32_MAX;
    for (size_t i = 0; i < hemesh->halfedge_count; i++) {
        if (hemesh->halfedges[i].vertex == vertex_idx) {
            he_start = i;
            break;
        }
    }

    if (he_start == UINT32_MAX) return max_length;

    /* Compute max normal variation */
    float max_variation = 0.0f;
    uint32_t he_current = he_start;
    do {
        uint32_t face = hemesh->halfedges[he_current].face;
        if (face != UINT32_MAX) {
            geomlib_vec3 face_normal = triangle_normal(mesh, face);
            float variation = 1.0f - geomlib_vec3_dot(normal, face_normal);
            if (variation > max_variation) max_variation = variation;
        }

        uint32_t twin = hemesh->halfedges[he_current].twin;
        if (twin == UINT32_MAX) break;
        he_current = hemesh->halfedges[twin].next;

        if (he_current == he_start) break;
    } while (he_current != UINT32_MAX);

    /* Map variation to edge length (high curvature = short edges) */
    float t = max_variation * 2.0f; /* Scale factor */
    t = (t < 0.0f) ? 0.0f : ((t > 1.0f) ? 1.0f : t);

    return min_length + (1.0f - t) * (max_length - min_length);
}

/* Default parameters */
geomlib_remeshing_params geomlib_remeshing_params_default(
    geomlib_remeshing_mode mode,
    float target_edge_length
) {
    geomlib_remeshing_params params;
    params.mode = mode;
    params.target_edge_length = target_edge_length;
    params.min_edge_length = target_edge_length * 0.5f;
    params.max_edge_length = target_edge_length * 2.0f;
    params.iterations = 10;
    params.feature_angle = 30.0f;
    params.preserve_boundaries = true;
    params.project_to_surface = false;
    params.split_threshold = 4.0f / 3.0f;
    params.collapse_threshold = 4.0f / 5.0f;
    params.smoothing_iterations = 1;
    return params;
}

/* Main remeshing function */
int geomlib_remesh(
    geomlib_trimesh *mesh,
    const geomlib_remeshing_params *params
) {
    if (!mesh || !params) return GEOMLIB_MESH_ERROR_INVALID;
    if (mesh->vertex_count == 0 || mesh->triangle_count == 0) return GEOMLIB_MESH_ERROR_INVALID;

    /* Build half-edge structure */
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Detect features */
    bool *is_feature = NULL;
    if (params->feature_angle > 0.0f) {
        is_feature = detect_features(mesh, &hemesh, params->feature_angle);
    }

    /* Store original vertex positions for projection */
    geomlib_vec3 *original_positions = NULL;
    if (params->project_to_surface) {
        original_positions = malloc(mesh->vertex_count * sizeof(geomlib_vec3));
        memcpy(original_positions, mesh->vertices, mesh->vertex_count * sizeof(geomlib_vec3));
    }

    /* Main remeshing loop */
    for (int iter = 0; iter < params->iterations; iter++) {
        /* 1. Edge Split */
        for (size_t face = 0; face < mesh->triangle_count; face++) {
            /* Check all three edges */
            for (int e = 0; e < 3; e++) {
                uint32_t v0 = mesh->indices[face * 3 + e];
                uint32_t v1 = mesh->indices[face * 3 + ((e + 1) % 3)];

                /* Skip feature edges */
                if (is_feature && is_feature[v0] && is_feature[v1]) continue;

                /* Compute target length */
                float target_len;
                if (params->mode == GEOMLIB_REMESH_ADAPTIVE) {
                    float len0 = compute_adaptive_edge_length(mesh, &hemesh, v0,
                                                               params->min_edge_length,
                                                               params->max_edge_length);
                    float len1 = compute_adaptive_edge_length(mesh, &hemesh, v1,
                                                               params->min_edge_length,
                                                               params->max_edge_length);
                    target_len = (len0 + len1) * 0.5f;
                } else {
                    target_len = params->target_edge_length;
                }

                float len = edge_length(mesh, v0, v1);
                if (len > target_len * params->split_threshold) {
                    split_edge(mesh, &hemesh, v0, v1);
                }
            }
        }

        /* 2. Edge Collapse */
        /* Note: We need to mark edges for collapse first to avoid modifying during iteration */
        bool collapsed_any = false;
        for (size_t he_idx = 0; he_idx < hemesh.halfedge_count; he_idx++) {
            uint32_t next = hemesh.halfedges[he_idx].next;
            if (next == UINT32_MAX) continue;

            uint32_t v0 = hemesh.halfedges[he_idx].vertex;
            uint32_t v1 = hemesh.halfedges[next].vertex;

            /* Skip feature edges and boundary edges */
            if (is_feature && is_feature[v0] && is_feature[v1]) continue;

            float target_len;
            if (params->mode == GEOMLIB_REMESH_ADAPTIVE) {
                float len0 = compute_adaptive_edge_length(mesh, &hemesh, v0,
                                                           params->min_edge_length,
                                                           params->max_edge_length);
                float len1 = compute_adaptive_edge_length(mesh, &hemesh, v1,
                                                           params->min_edge_length,
                                                           params->max_edge_length);
                target_len = (len0 + len1) * 0.5f;
            } else {
                target_len = params->target_edge_length;
            }

            float len = edge_length(mesh, v0, v1);
            if (len < target_len * params->collapse_threshold) {
                /* Compute new position (midpoint) */
                geomlib_vec3 p0 = mesh->vertices[v0];
                geomlib_vec3 p1 = mesh->vertices[v1];
                geomlib_vec3 new_pos = geomlib_vec3_mul(geomlib_vec3_add(p0, p1), 0.5f);

                /* Check validity */
                if (is_collapse_valid(mesh, &hemesh, v0, v1, new_pos)) {
                    /* Collapse by moving v0 to midpoint and removing v1 */
                    mesh->vertices[v0] = new_pos;

                    /* Update all references to v1 to point to v0 */
                    for (size_t f = 0; f < mesh->triangle_count; f++) {
                        for (int j = 0; j < 3; j++) {
                            if (mesh->indices[f * 3 + j] == v1) {
                                mesh->indices[f * 3 + j] = v0;
                            }
                        }
                    }

                    collapsed_any = true;
                }
            }
        }

        /* Remove degenerate triangles created by collapse */
        if (collapsed_any) {
            geomlib_remove_degenerate_triangles(mesh, 1e-10f);
            geomlib_halfedge_build(&hemesh, mesh);
        }

        /* 3. Tangential Smoothing */
        for (int smooth = 0; smooth < params->smoothing_iterations; smooth++) {
            tangential_smoothing(mesh, &hemesh, is_feature, params->preserve_boundaries);
        }

        /* 4. Project to surface (if enabled) */
        if (params->project_to_surface && original_positions) {
            /* Simple closest point projection (can be improved with BVH) */
            for (size_t v = 0; v < mesh->vertex_count; v++) {
                if (is_feature && is_feature[v]) continue;

                /* Find closest original vertex (very simple projection) */
                float min_dist = FLT_MAX;
                geomlib_vec3 closest = mesh->vertices[v];

                for (size_t ov = 0; ov < mesh->vertex_count && ov < 100; ov++) {
                    geomlib_vec3 diff = geomlib_vec3_sub(mesh->vertices[v], original_positions[ov]);
                    float dist = geomlib_vec3_length(diff);
                    if (dist < min_dist) {
                        min_dist = dist;
                        closest = original_positions[ov];
                    }
                }

                /* Blend with closest point */
                mesh->vertices[v] = geomlib_vec3_mul(
                    geomlib_vec3_add(mesh->vertices[v], closest), 0.5f
                );
            }
        }

        /* Rebuild half-edge for next iteration */
        if (iter < params->iterations - 1) {
            geomlib_halfedge_build(&hemesh, mesh);
        }
    }

    /* Cleanup */
    geomlib_halfedge_free(&hemesh);
    if (is_feature) free(is_feature);
    if (original_positions) free(original_positions);

    /* Recompute normals */
    geomlib_trimesh_compute_normals(mesh);

    return GEOMLIB_MESH_SUCCESS;
}

/* Simplified interfaces */
int geomlib_remesh_uniform(
    geomlib_trimesh *mesh,
    float target_edge_length,
    int iterations
) {
    geomlib_remeshing_params params = geomlib_remeshing_params_default(
        GEOMLIB_REMESH_UNIFORM, target_edge_length
    );
    params.iterations = iterations;
    return geomlib_remesh(mesh, &params);
}

int geomlib_remesh_adaptive(
    geomlib_trimesh *mesh,
    float min_edge_length,
    float max_edge_length,
    int iterations
) {
    geomlib_remeshing_params params = geomlib_remeshing_params_default(
        GEOMLIB_REMESH_ADAPTIVE, (min_edge_length + max_edge_length) * 0.5f
    );
    params.min_edge_length = min_edge_length;
    params.max_edge_length = max_edge_length;
    params.iterations = iterations;
    return geomlib_remesh(mesh, &params);
}

/* Compute target edge length */
float geomlib_compute_target_edge_length(
    const geomlib_trimesh *mesh,
    size_t target_vertex_count
) {
    if (!mesh || target_vertex_count == 0) return 0.0f;

    /* Compute surface area */
    float area = 0.0f;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t i0 = mesh->indices[i * 3 + 0];
        uint32_t i1 = mesh->indices[i * 3 + 1];
        uint32_t i2 = mesh->indices[i * 3 + 2];

        geomlib_vec3 p0 = mesh->vertices[i0];
        geomlib_vec3 p1 = mesh->vertices[i1];
        geomlib_vec3 p2 = mesh->vertices[i2];

        geomlib_vec3 e1 = geomlib_vec3_sub(p1, p0);
        geomlib_vec3 e2 = geomlib_vec3_sub(p2, p0);
        geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);
        area += geomlib_vec3_length(cross) * 0.5f;
    }

    /* Estimate edge length from area and target vertex count */
    /* For equilateral triangles: A = sqrt(3)/4 * L^2 */
    /* For mesh: total_area ≈ num_triangles * sqrt(3)/4 * L^2 */
    /* Euler: V - E + F = 2, E ≈ 3V, F ≈ 2V for large meshes */
    float num_triangles = target_vertex_count * 2.0f;
    float L_squared = area / (num_triangles * 0.433f); /* sqrt(3)/4 ≈ 0.433 */

    return sqrtf(L_squared);
}

/* Compute mean edge length */
float geomlib_compute_mean_edge_length(const geomlib_trimesh *mesh) {
    if (!mesh || mesh->triangle_count == 0) return 0.0f;

    float total_length = 0.0f;
    int edge_count = 0;

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t i0 = mesh->indices[i * 3 + 0];
        uint32_t i1 = mesh->indices[i * 3 + 1];
        uint32_t i2 = mesh->indices[i * 3 + 2];

        total_length += edge_length(mesh, i0, i1);
        total_length += edge_length(mesh, i1, i2);
        total_length += edge_length(mesh, i2, i0);
        edge_count += 3;
    }

    return total_length / edge_count;
}

/* Compute mesh quality statistics */
void geomlib_compute_remeshing_stats(
    const geomlib_trimesh *mesh,
    geomlib_remeshing_stats *stats
) {
    if (!mesh || !stats) return;

    memset(stats, 0, sizeof(geomlib_remeshing_stats));
    stats->min_edge_length = FLT_MAX;
    stats->min_angle = 180.0f;
    stats->max_angle = 0.0f;

    float total_length = 0.0f;
    int edge_count = 0;

    /* Analyze triangles */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t i0 = mesh->indices[i * 3 + 0];
        uint32_t i1 = mesh->indices[i * 3 + 1];
        uint32_t i2 = mesh->indices[i * 3 + 2];

        /* Edge lengths */
        float l01 = edge_length(mesh, i0, i1);
        float l12 = edge_length(mesh, i1, i2);
        float l20 = edge_length(mesh, i2, i0);

        if (l01 < stats->min_edge_length) stats->min_edge_length = l01;
        if (l12 < stats->min_edge_length) stats->min_edge_length = l12;
        if (l20 < stats->min_edge_length) stats->min_edge_length = l20;

        if (l01 > stats->max_edge_length) stats->max_edge_length = l01;
        if (l12 > stats->max_edge_length) stats->max_edge_length = l12;
        if (l20 > stats->max_edge_length) stats->max_edge_length = l20;

        total_length += l01 + l12 + l20;
        edge_count += 3;

        /* Triangle angles using law of cosines */
        float angle0 = acosf((l01*l01 + l20*l20 - l12*l12) / (2.0f * l01 * l20)) * 180.0f / M_PI;
        float angle1 = acosf((l01*l01 + l12*l12 - l20*l20) / (2.0f * l01 * l12)) * 180.0f / M_PI;
        float angle2 = acosf((l12*l12 + l20*l20 - l01*l01) / (2.0f * l12 * l20)) * 180.0f / M_PI;

        if (angle0 < stats->min_angle) stats->min_angle = angle0;
        if (angle1 < stats->min_angle) stats->min_angle = angle1;
        if (angle2 < stats->min_angle) stats->min_angle = angle2;

        if (angle0 > stats->max_angle) stats->max_angle = angle0;
        if (angle1 > stats->max_angle) stats->max_angle = angle1;
        if (angle2 > stats->max_angle) stats->max_angle = angle2;
    }

    stats->mean_edge_length = total_length / edge_count;

    /* Count irregular vertices */
    geomlib_halfedge_mesh hemesh;
    if (geomlib_halfedge_build(&hemesh, mesh) == GEOMLIB_MESH_SUCCESS) {
        for (size_t v = 0; v < mesh->vertex_count; v++) {
            int valence = get_vertex_valence(&hemesh, v);
            if (valence != 6) stats->irregular_vertices++;
        }
        geomlib_halfedge_free(&hemesh);
    }
}
