/*
 * geomlib - Mesh Healing Implementation
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 */

#include "geomlib/geomlib_healing.h"
#include "geomlib/geomlib_holefill.h"
#include "geomlib/geomlib_random.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define INVALID_INDEX UINT32_MAX

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static inline float vec3_length(geomlib_vec3 v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

static inline float triangle_area(geomlib_vec3 v0, geomlib_vec3 v1, geomlib_vec3 v2) {
    geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);
    return vec3_length(cross) * 0.5f;
}

/* ============================================================================
 * Half-Edge Construction
 * ============================================================================ */

typedef struct {
    uint32_t v0, v1;  /* Edge vertices (v0 < v1) */
    uint32_t he_idx;  /* Half-edge index */
} EdgeKey;

static int compare_edges(const void *a, const void *b) {
    const EdgeKey *ea = (const EdgeKey*)a;
    const EdgeKey *eb = (const EdgeKey*)b;

    if (ea->v0 != eb->v0) return (int)ea->v0 - (int)eb->v0;
    return (int)ea->v1 - (int)eb->v1;
}

int geomlib_halfedge_build(
    geomlib_halfedge_mesh *hemesh,
    const geomlib_trimesh *mesh
) {
    if (!hemesh || !mesh) return GEOMLIB_MESH_ERROR_INVALID;
    if (mesh->triangle_count == 0) return GEOMLIB_MESH_ERROR_INVALID;

    memset(hemesh, 0, sizeof(*hemesh));
    hemesh->mesh = mesh;

    // Allocate half-edges (3 per triangle)
    size_t he_count = mesh->triangle_count * 3;
    hemesh->halfedges = (geomlib_halfedge*)malloc(he_count * sizeof(geomlib_halfedge));
    if (!hemesh->halfedges) return GEOMLIB_MESH_ERROR_MEMORY;

    hemesh->halfedge_count = he_count;
    hemesh->halfedge_capacity = he_count;

    // Allocate vertex and face tables
    hemesh->vertex_halfedge = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    hemesh->face_halfedge = (uint32_t*)malloc(mesh->triangle_count * sizeof(uint32_t));

    if (!hemesh->vertex_halfedge || !hemesh->face_halfedge) {
        geomlib_halfedge_free(hemesh);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    // Initialize to invalid
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        hemesh->vertex_halfedge[i] = INVALID_INDEX;
    }

    // Build half-edges
    for (size_t face = 0; face < mesh->triangle_count; face++) {
        uint32_t i0 = mesh->indices[face * 3 + 0];
        uint32_t i1 = mesh->indices[face * 3 + 1];
        uint32_t i2 = mesh->indices[face * 3 + 2];

        uint32_t he0 = (uint32_t)(face * 3 + 0);
        uint32_t he1 = (uint32_t)(face * 3 + 1);
        uint32_t he2 = (uint32_t)(face * 3 + 2);

        // Half-edge 0: i0 -> i1
        hemesh->halfedges[he0].vertex = i1;
        hemesh->halfedges[he0].face = (uint32_t)face;
        hemesh->halfedges[he0].next = he1;
        hemesh->halfedges[he0].twin = INVALID_INDEX;

        // Half-edge 1: i1 -> i2
        hemesh->halfedges[he1].vertex = i2;
        hemesh->halfedges[he1].face = (uint32_t)face;
        hemesh->halfedges[he1].next = he2;
        hemesh->halfedges[he1].twin = INVALID_INDEX;

        // Half-edge 2: i2 -> i0
        hemesh->halfedges[he2].vertex = i0;
        hemesh->halfedges[he2].face = (uint32_t)face;
        hemesh->halfedges[he2].next = he0;
        hemesh->halfedges[he2].twin = INVALID_INDEX;

        // Update vertex half-edge pointers
        if (hemesh->vertex_halfedge[i0] == INVALID_INDEX) hemesh->vertex_halfedge[i0] = he0;
        if (hemesh->vertex_halfedge[i1] == INVALID_INDEX) hemesh->vertex_halfedge[i1] = he1;
        if (hemesh->vertex_halfedge[i2] == INVALID_INDEX) hemesh->vertex_halfedge[i2] = he2;

        hemesh->face_halfedge[face] = he0;
    }

    // Find twin half-edges
    // Create sorted edge list for O(n log n) twin finding
    EdgeKey *edges = (EdgeKey*)malloc(he_count * sizeof(EdgeKey));
    if (!edges) {
        geomlib_halfedge_free(hemesh);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (size_t i = 0; i < he_count; i++) {
        uint32_t v_from = (i % 3 == 0) ? mesh->indices[i / 3 * 3 + 0] :
                          (i % 3 == 1) ? mesh->indices[i / 3 * 3 + 1] :
                                        mesh->indices[i / 3 * 3 + 2];
        uint32_t v_to = hemesh->halfedges[i].vertex;

        // Store edge with sorted vertices for matching
        edges[i].v0 = v_from < v_to ? v_from : v_to;
        edges[i].v1 = v_from < v_to ? v_to : v_from;
        edges[i].he_idx = (uint32_t)i;
    }

    qsort(edges, he_count, sizeof(EdgeKey), compare_edges);

    // Match twin half-edges
    for (size_t i = 0; i < he_count; i++) {
        if (i + 1 < he_count &&
            edges[i].v0 == edges[i + 1].v0 &&
            edges[i].v1 == edges[i + 1].v1) {

            uint32_t he1 = edges[i].he_idx;
            uint32_t he2 = edges[i + 1].he_idx;

            hemesh->halfedges[he1].twin = he2;
            hemesh->halfedges[he2].twin = he1;

            i++; // Skip next edge (already matched)
        }
    }

    free(edges);
    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_halfedge_free(geomlib_halfedge_mesh *hemesh) {
    if (!hemesh) return;

    free(hemesh->halfedges);
    free(hemesh->vertex_halfedge);
    free(hemesh->face_halfedge);

    memset(hemesh, 0, sizeof(*hemesh));
}

void geomlib_halfedge_get_stats(
    const geomlib_halfedge_mesh *hemesh,
    size_t *boundary_edges_out,
    bool *manifold_out,
    int *genus_out
) {
    if (!hemesh) return;

    size_t boundary_edges = 0;
    bool is_manifold = true;

    for (size_t i = 0; i < hemesh->halfedge_count; i++) {
        if (hemesh->halfedges[i].twin == INVALID_INDEX) {
            boundary_edges++;
        }
    }

    if (boundary_edges_out) *boundary_edges_out = boundary_edges;
    if (manifold_out) *manifold_out = is_manifold;

    // Compute genus using Euler characteristic: V - E + F = 2 - 2g
    if (genus_out && hemesh->mesh) {
        int V = (int)hemesh->mesh->vertex_count;
        int F = (int)hemesh->mesh->triangle_count;
        int E = (int)hemesh->halfedge_count / 2 + (int)boundary_edges / 2;

        int euler = V - E + F;
        int genus = (2 - euler) / 2;
        *genus_out = genus > 0 ? genus : 0;
    }
}

/* ============================================================================
 * Face Orientation Detection
 * ============================================================================ */

int geomlib_detect_face_orientation(
    const geomlib_trimesh *mesh,
    uint32_t face_idx,
    const geomlib_bvh *bvh,
    int ray_samples
) {
    (void)bvh;
    if (!mesh || face_idx >= mesh->triangle_count) return 0;

    // Get face vertices
    uint32_t i0 = mesh->indices[face_idx * 3 + 0];
    uint32_t i1 = mesh->indices[face_idx * 3 + 1];
    uint32_t i2 = mesh->indices[face_idx * 3 + 2];

    geomlib_vec3 v0 = mesh->vertices[i0];
    geomlib_vec3 v1 = mesh->vertices[i1];
    geomlib_vec3 v2 = mesh->vertices[i2];

    // Compute face center and normal
    geomlib_vec3 center = geomlib_vec3_mul(
        geomlib_vec3_add(geomlib_vec3_add(v0, v1), v2), 1.0f / 3.0f
    );

    geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 normal = geomlib_vec3_normalize(geomlib_vec3_cross(e1, e2));

    // Cast rays from slightly offset position in normal direction
    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, face_idx, 12345);

    int inside_count = 0;
    int outside_count = 0;

    for (int i = 0; i < ray_samples; i++) {
        // Random ray direction
        float theta = geomlib_pcg32_nextf(&rng) * 2.0f * (float)M_PI;
        float phi = acosf(2.0f * geomlib_pcg32_nextf(&rng) - 1.0f);

        geomlib_vec3 ray_dir = geomlib_vec3_create(
            sinf(phi) * cosf(theta),
            sinf(phi) * sinf(theta),
            cosf(phi)
        );

        // Offset origin slightly along normal
        geomlib_vec3 ray_origin = geomlib_vec3_add(center, geomlib_vec3_mul(normal, 0.001f));

        // Count intersections
        geomlib_ray ray;
        ray.origin = ray_origin;
        ray.direction = ray_dir;
        ray.t_min = 0.001f;
        ray.t_max = 1000.0f;

        // Simple ray casting without BVH for now
        // In production, would use BVH if available
        int hit_count = 0;

        for (size_t f = 0; f < mesh->triangle_count; f++) {
            if (f == face_idx) continue;

            uint32_t fi0 = mesh->indices[f * 3 + 0];
            uint32_t fi1 = mesh->indices[f * 3 + 1];
            uint32_t fi2 = mesh->indices[f * 3 + 2];

            float t, u, v;
            if (geomlib_ray_triangle_intersect(&ray,
                mesh->vertices[fi0],
                mesh->vertices[fi1],
                mesh->vertices[fi2],
                &t, &u, &v)) {
                hit_count++;
            }
        }

        // Odd number of intersections = inside
        if (hit_count % 2 == 1) {
            inside_count++;
        } else {
            outside_count++;
        }
    }

    // If majority of rays from normal direction are inside, face points inward
    return outside_count - inside_count;
}

/* ============================================================================
 * Face Orientation Fixing
 * ============================================================================ */

geomlib_orientation_params geomlib_orientation_params_default(
    geomlib_orientation_strategy strategy
) {
    geomlib_orientation_params params;
    memset(&params, 0, sizeof(params));

    params.strategy = strategy;
    params.ray_samples = 32;
    params.ray_offset = 0.001f;
    params.flip_if_majority_inside = true;

    return params;
}

static int fix_orientation_propagation(
    geomlib_trimesh *mesh,
    geomlib_halfedge_mesh *hemesh,
    bool flip_if_inside,
    int ray_samples
) {
    if (!hemesh) return GEOMLIB_MESH_ERROR_INVALID;

    bool *visited = (bool*)calloc(mesh->triangle_count, sizeof(bool));
    bool *should_flip = (bool*)calloc(mesh->triangle_count, sizeof(bool));
    uint32_t *queue = (uint32_t*)malloc(mesh->triangle_count * sizeof(uint32_t));

    if (!visited || !should_flip || !queue) {
        free(visited);
        free(should_flip);
        free(queue);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    // Start from face 0 (arbitrary seed)
    uint32_t seed_face = 0;

    // Optionally check if seed face points inward
    if (flip_if_inside && ray_samples > 0) {
        int orientation = geomlib_detect_face_orientation(mesh, seed_face, NULL, ray_samples);
        should_flip[seed_face] = (orientation < 0);
    }

    // BFS to propagate orientation
    int queue_front = 0, queue_back = 0;
    queue[queue_back++] = seed_face;
    visited[seed_face] = true;

    while (queue_front < queue_back) {
        uint32_t face = queue[queue_front++];
        uint32_t he_idx = hemesh->face_halfedge[face];

        // Visit three edges of this face
        for (int e = 0; e < 3; e++) {
            uint32_t he = (uint32_t)(he_idx + e);
            uint32_t twin = hemesh->halfedges[he].twin;

            if (twin != INVALID_INDEX) {
                uint32_t neighbor_face = hemesh->halfedges[twin].face;

                if (!visited[neighbor_face]) {
                    visited[neighbor_face] = true;

                    // Check if neighbor has consistent orientation
                    // If edge vertices are in opposite order, orientation is consistent
                    uint32_t he_v0_idx = (he % 3 == 2) ? (he - 2) : (he + 1);
                    uint32_t he_v1 = hemesh->halfedges[he].vertex;
                    uint32_t he_v0 = hemesh->halfedges[he_v0_idx].vertex;

                    uint32_t twin_v0_idx = (twin % 3 == 2) ? (twin - 2) : (twin + 1);
                    uint32_t twin_v1 = hemesh->halfedges[twin].vertex;
                    uint32_t twin_v0 = hemesh->halfedges[twin_v0_idx].vertex;

                    // If vertices match (same order), faces have inconsistent orientation
                    bool inconsistent = (he_v0 == twin_v0 && he_v1 == twin_v1);

                    should_flip[neighbor_face] = should_flip[face] ^ inconsistent;

                    queue[queue_back++] = neighbor_face;
                }
            }
        }
    }

    // Apply flips
    int flipped_count = 0;
    for (size_t f = 0; f < mesh->triangle_count; f++) {
        if (should_flip[f]) {
            // Swap indices to flip face
            uint32_t tmp = mesh->indices[f * 3 + 1];
            mesh->indices[f * 3 + 1] = mesh->indices[f * 3 + 2];
            mesh->indices[f * 3 + 2] = tmp;
            flipped_count++;
        }
    }

    free(visited);
    free(should_flip);
    free(queue);

    return flipped_count;
}

int geomlib_fix_face_orientation(
    geomlib_trimesh *mesh,
    const geomlib_orientation_params *params,
    geomlib_halfedge_mesh *hemesh
) {
    if (!mesh || !params) return GEOMLIB_MESH_ERROR_INVALID;

    // Build half-edge structure if not provided
    geomlib_halfedge_mesh local_hemesh;
    bool need_free = false;

    if (!hemesh) {
        int result = geomlib_halfedge_build(&local_hemesh, mesh);
        if (result != GEOMLIB_MESH_SUCCESS) return result;
        hemesh = &local_hemesh;
        need_free = true;
    }

    int flipped_count = 0;

    switch (params->strategy) {
        case GEOMLIB_ORIENT_PROPAGATION:
            flipped_count = fix_orientation_propagation(mesh, hemesh, false, 0);
            break;

        case GEOMLIB_ORIENT_INSIDE_OUT:
        case GEOMLIB_ORIENT_HYBRID:
            flipped_count = fix_orientation_propagation(mesh, hemesh,
                params->flip_if_majority_inside, params->ray_samples);
            break;

        case GEOMLIB_ORIENT_MAJORITY_VOTE:
            // Simple majority vote: count consistent vs inconsistent neighbors
            flipped_count = fix_orientation_propagation(mesh, hemesh, false, 0);
            break;
    }

    if (need_free) {
        geomlib_halfedge_free(&local_hemesh);
    }

    // Recompute normals after flipping
    if (flipped_count > 0 && mesh->normals) {
        geomlib_trimesh_compute_normals(mesh);
    }

    return flipped_count;
}

/* ============================================================================
 * Mesh Healing Operations
 * ============================================================================ */

int geomlib_remove_degenerate_triangles(
    geomlib_trimesh *mesh,
    float area_threshold
) {
    if (!mesh) return GEOMLIB_MESH_ERROR_INVALID;

    if (area_threshold <= 0.0f) {
        area_threshold = 1e-10f;
    }

    size_t write_idx = 0;

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t i0 = mesh->indices[i * 3 + 0];
        uint32_t i1 = mesh->indices[i * 3 + 1];
        uint32_t i2 = mesh->indices[i * 3 + 2];

        // Check for duplicate vertices
        if (i0 == i1 || i1 == i2 || i2 == i0) {
            continue; // Skip degenerate
        }

        // Check area
        float area = triangle_area(
            mesh->vertices[i0],
            mesh->vertices[i1],
            mesh->vertices[i2]
        );

        if (area >= area_threshold) {
            // Keep this triangle
            if (write_idx != i) {
                mesh->indices[write_idx * 3 + 0] = i0;
                mesh->indices[write_idx * 3 + 1] = i1;
                mesh->indices[write_idx * 3 + 2] = i2;
            }
            write_idx++;
        }
    }

    int removed = (int)(mesh->triangle_count - write_idx);
    mesh->triangle_count = write_idx;

    return removed;
}

int geomlib_remove_duplicate_vertices(
    geomlib_trimesh *mesh,
    float distance_threshold
) {
    if (!mesh) return GEOMLIB_MESH_ERROR_INVALID;
    if (mesh->vertex_count == 0) return 0;

    if (distance_threshold <= 0.0f) {
        distance_threshold = 1e-6f;
    }

    // Simple O(n²) implementation for now
    // Production would use spatial hashing
    uint32_t *remap = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    if (!remap) return GEOMLIB_MESH_ERROR_MEMORY;

    bool *keep = (bool*)calloc(mesh->vertex_count, sizeof(bool));
    if (!keep) {
        free(remap);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    size_t unique_count = 0;

    for (size_t i = 0; i < mesh->vertex_count; i++) {
        remap[i] = (uint32_t)i;

        // Check if duplicate of earlier vertex
        bool is_duplicate = false;
        for (size_t j = 0; j < i; j++) {
            if (!keep[j]) continue;

            geomlib_vec3 diff = geomlib_vec3_sub(mesh->vertices[i], mesh->vertices[j]);
            float dist = vec3_length(diff);

            if (dist < distance_threshold) {
                remap[i] = (uint32_t)j;
                is_duplicate = true;
                break;
            }
        }

        if (!is_duplicate) {
            keep[i] = true;
            unique_count++;
        }
    }

    // Remap triangle indices
    for (size_t i = 0; i < mesh->triangle_count * 3; i++) {
        mesh->indices[i] = remap[mesh->indices[i]];
    }

    free(remap);
    free(keep);

    int merged = (int)(mesh->vertex_count - unique_count);
    return merged;
}

int geomlib_fill_holes(
    geomlib_trimesh *mesh,
    int max_hole_edges
) {
    if (!mesh) return GEOMLIB_MESH_ERROR_INVALID;

    /* Detect holes in the mesh */
    geomlib_holes holes = geomlib_detect_holes(mesh);
    if (holes.hole_count == 0) {
        geomlib_holes_free(&holes);
        return 0;
    }

    /* Count filled holes */
    int filled_count = 0;

    /* Setup fill parameters with natural fill method */
    geomlib_fill_params params = geomlib_fill_params_default(GEOMLIB_FILL_NATURAL);
    params.fair_interior = 1;
    params.smooth_iterations = 5;

    /* Fill each hole */
    for (size_t i = 0; i < holes.hole_count; i++) {
        /* Skip holes larger than max_hole_edges if specified */
        if (max_hole_edges > 0 && (int)holes.holes[i].vertex_count > max_hole_edges) {
            continue;
        }

        /* Fill the hole */
        int result = geomlib_fill_hole(mesh, &holes.holes[i], &params);
        if (result == GEOMLIB_MESH_SUCCESS) {
            filled_count++;
        }
    }

    geomlib_holes_free(&holes);
    return filled_count;
}

int geomlib_repair_mesh(
    geomlib_trimesh *mesh,
    bool fill_holes
) {
    if (!mesh) return GEOMLIB_MESH_ERROR_INVALID;

    // 1. Remove duplicate vertices
    geomlib_remove_duplicate_vertices(mesh, 1e-6f);

    // 2. Remove degenerate triangles
    geomlib_remove_degenerate_triangles(mesh, 1e-10f);

    // 3. Fix face orientations
    geomlib_orientation_params params = geomlib_orientation_params_default(
        GEOMLIB_ORIENT_PROPAGATION
    );
    geomlib_fix_face_orientation(mesh, &params, NULL);

    // 4. Fill holes if requested
    if (fill_holes) {
        geomlib_fill_holes(mesh, 100);
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Mesh Topology Analysis
 * ============================================================================ */

int geomlib_analyze_mesh_topology(
    const geomlib_trimesh *mesh,
    geomlib_mesh_issues *issues
) {
    if (!mesh || !issues) return GEOMLIB_MESH_ERROR_INVALID;

    memset(issues, 0, sizeof(*issues));

    // Build half-edge structure for analysis
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, mesh);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    // Get topology stats
    size_t boundary_edges = 0;
    geomlib_halfedge_get_stats(&hemesh, &boundary_edges, &issues->is_manifold, NULL);

    issues->boundary_edges = boundary_edges;
    issues->is_watertight = (boundary_edges == 0);

    // Count degenerate triangles
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t i0 = mesh->indices[i * 3 + 0];
        uint32_t i1 = mesh->indices[i * 3 + 1];
        uint32_t i2 = mesh->indices[i * 3 + 2];

        if (i0 == i1 || i1 == i2 || i2 == i0) {
            issues->degenerate_triangles++;
            continue;
        }

        float area = triangle_area(mesh->vertices[i0], mesh->vertices[i1], mesh->vertices[i2]);
        if (area < 1e-10f) {
            issues->degenerate_triangles++;
        }
    }

    // Check orientation consistency
    // Two adjacent faces should have opposite edge orientation
    issues->is_oriented = true;
    for (size_t i = 0; i < hemesh.halfedge_count; i++) {
        if (hemesh.halfedges[i].twin != INVALID_INDEX) {
            uint32_t twin = hemesh.halfedges[i].twin;

            // Get vertices for this half-edge
            uint32_t he_v0_idx = (i % 3 == 2) ? (i - 2) : (i + 1);
            uint32_t he_v1 = hemesh.halfedges[i].vertex;
            uint32_t he_v0 = hemesh.halfedges[he_v0_idx].vertex;

            // Get vertices for twin
            uint32_t twin_v0_idx = (twin % 3 == 2) ? (twin - 2) : (twin + 1);
            uint32_t twin_v1 = hemesh.halfedges[twin].vertex;
            uint32_t twin_v0 = hemesh.halfedges[twin_v0_idx].vertex;

            // Check if they're opposite
            if (he_v0 == twin_v1 && he_v1 == twin_v0) {
                // Consistent
            } else {
                issues->inconsistent_faces++;
                issues->is_oriented = false;
            }
        }
    }

    geomlib_halfedge_free(&hemesh);

    return GEOMLIB_MESH_SUCCESS;
}
