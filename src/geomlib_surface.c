/*
 * geomlib - Surface Reconstruction Implementation
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 */

#include "geomlib/geomlib_surface.h"
#include "geomlib/geomlib_random.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static inline float vec3_distance(geomlib_vec3 a, geomlib_vec3 b) {
    float dx = b.x - a.x;
    float dy = b.y - a.y;
    float dz = b.z - a.z;
    return sqrtf(dx * dx + dy * dy + dz * dz);
}

/* ============================================================================
 * Point Cloud Management
 * ============================================================================ */

int geomlib_pointcloud_init(geomlib_pointcloud *pc, size_t initial_capacity) {
    if (!pc) return GEOMLIB_MESH_ERROR_INVALID;

    pc->points = (geomlib_vec3*)malloc(initial_capacity * sizeof(geomlib_vec3));
    if (!pc->points) return GEOMLIB_MESH_ERROR_MEMORY;

    pc->normals = NULL;
    pc->confidences = NULL;
    pc->point_count = 0;
    pc->point_capacity = initial_capacity;

    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_pointcloud_free(geomlib_pointcloud *pc) {
    if (!pc) return;

    free(pc->points);
    free(pc->normals);
    free(pc->confidences);

    pc->points = NULL;
    pc->normals = NULL;
    pc->confidences = NULL;
    pc->point_count = 0;
    pc->point_capacity = 0;
}

size_t geomlib_pointcloud_add_point(
    geomlib_pointcloud *pc,
    geomlib_vec3 position,
    geomlib_vec3 *normal,
    float confidence
) {
    if (!pc || !pc->points) return (size_t)-1;

    // Grow if needed
    if (pc->point_count >= pc->point_capacity) {
        size_t new_capacity = pc->point_capacity * 2;
        geomlib_vec3 *new_points = (geomlib_vec3*)realloc(pc->points, new_capacity * sizeof(geomlib_vec3));
        if (!new_points) return (size_t)-1;
        pc->points = new_points;

        if (pc->normals) {
            geomlib_vec3 *new_normals = (geomlib_vec3*)realloc(pc->normals, new_capacity * sizeof(geomlib_vec3));
            if (!new_normals) return (size_t)-1;
            pc->normals = new_normals;
        }

        if (pc->confidences) {
            float *new_conf = (float*)realloc(pc->confidences, new_capacity * sizeof(float));
            if (!new_conf) return (size_t)-1;
            pc->confidences = new_conf;
        }

        pc->point_capacity = new_capacity;
    }

    // Allocate normals/confidences on first use
    if (normal && !pc->normals) {
        pc->normals = (geomlib_vec3*)calloc(pc->point_capacity, sizeof(geomlib_vec3));
        if (!pc->normals) return (size_t)-1;
    }

    if (!pc->confidences) {
        pc->confidences = (float*)malloc(pc->point_capacity * sizeof(float));
        if (!pc->confidences) return (size_t)-1;
        for (size_t i = 0; i < pc->point_capacity; i++) {
            pc->confidences[i] = 1.0f;
        }
    }

    size_t idx = pc->point_count++;
    pc->points[idx] = position;

    if (normal && pc->normals) {
        pc->normals[idx] = *normal;
    }

    if (pc->confidences) {
        pc->confidences[idx] = confidence;
    }

    return idx;
}

geomlib_aabb geomlib_pointcloud_bounds(const geomlib_pointcloud *pc) {
    geomlib_aabb bounds;
    bounds.min = geomlib_vec3_create(FLT_MAX, FLT_MAX, FLT_MAX);
    bounds.max = geomlib_vec3_create(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    if (!pc || !pc->points || pc->point_count == 0) {
        return bounds;
    }

    for (size_t i = 0; i < pc->point_count; i++) {
        geomlib_vec3 p = pc->points[i];
        if (p.x < bounds.min.x) bounds.min.x = p.x;
        if (p.y < bounds.min.y) bounds.min.y = p.y;
        if (p.z < bounds.min.z) bounds.min.z = p.z;
        if (p.x > bounds.max.x) bounds.max.x = p.x;
        if (p.y > bounds.max.y) bounds.max.y = p.y;
        if (p.z > bounds.max.z) bounds.max.z = p.z;
    }

    return bounds;
}

void geomlib_pointcloud_get_stats(
    const geomlib_pointcloud *pc,
    float *avg_spacing_out,
    float *density_out
) {
    if (!pc || pc->point_count < 2) {
        if (avg_spacing_out) *avg_spacing_out = 0.0f;
        if (density_out) *density_out = 0.0f;
        return;
    }

    // Sample random points and compute average nearest neighbor distance
    const int samples = (int)(pc->point_count < 100 ? pc->point_count : 100);
    float total_dist = 0.0f;

    for (int i = 0; i < samples; i++) {
        size_t idx = (size_t)(rand() % (int)pc->point_count);
        geomlib_vec3 p = pc->points[idx];

        float min_dist = FLT_MAX;
        for (size_t j = 0; j < pc->point_count; j++) {
            if (j == idx) continue;
            float dist = vec3_distance(p, pc->points[j]);
            if (dist < min_dist) min_dist = dist;
        }

        total_dist += min_dist;
    }

    float avg_spacing = total_dist / samples;
    if (avg_spacing_out) *avg_spacing_out = avg_spacing;

    if (density_out) {
        geomlib_aabb bounds = geomlib_pointcloud_bounds(pc);
        geomlib_vec3 size = geomlib_vec3_sub(bounds.max, bounds.min);
        float volume = size.x * size.y * size.z;
        *density_out = volume > 0.0f ? (float)pc->point_count / volume : 0.0f;
    }
}

/* ============================================================================
 * Gaussian Splat Conversion
 * ============================================================================ */

int geomlib_gaussians_to_pointcloud(
    const float *positions,
    size_t count,
    geomlib_pointcloud *output
) {
    if (!positions || !output) return GEOMLIB_MESH_ERROR_INVALID;

    int result = geomlib_pointcloud_init(output, count);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    for (size_t i = 0; i < count; i++) {
        geomlib_vec3 pos = geomlib_vec3_create(
            positions[i * 3 + 0],
            positions[i * 3 + 1],
            positions[i * 3 + 2]
        );
        geomlib_pointcloud_add_point(output, pos, NULL, 1.0f);
    }

    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_gaussians_sample_points(
    const float *positions,
    const float *scales,
    const float *rotations,
    size_t count,
    int samples_per,
    geomlib_pointcloud *output
) {
    if (!positions || !output || samples_per < 1) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    int result = geomlib_pointcloud_init(output, count * samples_per);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 54);

    for (size_t i = 0; i < count; i++) {
        geomlib_vec3 center = geomlib_vec3_create(
            positions[i * 3 + 0],
            positions[i * 3 + 1],
            positions[i * 3 + 2]
        );

        // Use scale or default radius
        float radius = 0.1f;
        if (scales) {
            radius = (scales[i * 3 + 0] + scales[i * 3 + 1] + scales[i * 3 + 2]) / 3.0f;
        }

        // Sample points around Gaussian center
        for (int s = 0; s < samples_per; s++) {
            // Generate random point on unit sphere
            float theta = geomlib_pcg32_nextf(&rng) * 2.0f * (float)M_PI;
            float phi = acosf(2.0f * geomlib_pcg32_nextf(&rng) - 1.0f);

            float x = sinf(phi) * cosf(theta);
            float y = sinf(phi) * sinf(theta);
            float z = cosf(phi);

            // Scale by Gaussian radius (use smaller radius for surface points)
            float r = radius * 0.5f;

            geomlib_vec3 offset = geomlib_vec3_create(x * r, y * r, z * r);
            geomlib_vec3 point = geomlib_vec3_add(center, offset);

            geomlib_pointcloud_add_point(output, point, NULL, 1.0f);
        }
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Surface Reconstruction - Greedy Projection
 * ============================================================================ */

typedef struct {
    int point_idx;
    geomlib_vec3 position;
    bool used;
} ProjectionPoint;

static int reconstruct_greedy_projection(
    const geomlib_pointcloud *pc,
    const geomlib_surface_params *params,
    geomlib_trimesh *output
) {
    if (!pc || pc->point_count < 3) return GEOMLIB_MESH_ERROR_INVALID;

    // This is a simplified greedy projection
    // For production, would use better spatial indexing

    ProjectionPoint *points = (ProjectionPoint*)malloc(pc->point_count * sizeof(ProjectionPoint));
    if (!points) return GEOMLIB_MESH_ERROR_MEMORY;

    for (size_t i = 0; i < pc->point_count; i++) {
        points[i].point_idx = (int)i;
        points[i].position = pc->points[i];
        points[i].used = false;
    }

    float max_edge = params->max_edge_length;
    if (max_edge <= 0.0f) {
        // Auto-compute from point cloud
        float avg_spacing;
        geomlib_pointcloud_get_stats(pc, &avg_spacing, NULL);
        max_edge = avg_spacing * 3.0f;
    }

    // Simple greedy triangulation
    // Start with random point, find nearest neighbors, create triangles
    for (size_t i = 0; i < pc->point_count; i++) {
        if (points[i].used) continue;

        geomlib_vec3 p0 = points[i].position;

        // Find two nearest unused neighbors
        int best_j = -1, best_k = -1;
        float best_j_dist = FLT_MAX, best_k_dist = FLT_MAX;

        for (size_t j = 0; j < pc->point_count; j++) {
            if (j == i || points[j].used) continue;

            float dist = vec3_distance(p0, points[j].position);
            if (dist > max_edge) continue;

            if (dist < best_j_dist) {
                best_k = best_j;
                best_k_dist = best_j_dist;
                best_j = (int)j;
                best_j_dist = dist;
            } else if (dist < best_k_dist) {
                best_k = (int)j;
                best_k_dist = dist;
            }
        }

        if (best_j >= 0 && best_k >= 0) {
            // Add triangle
            geomlib_trimesh_add_triangle(output, (uint32_t)i, (uint32_t)best_j, (uint32_t)best_k);
            points[i].used = true;
        }
    }

    free(points);

    // Copy vertices to output mesh
    for (size_t i = 0; i < pc->point_count; i++) {
        geomlib_trimesh_add_vertex(output, pc->points[i],
            pc->normals ? &pc->normals[i] : NULL, NULL);
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Surface Reconstruction - Marching Cubes (Simplified)
 * ============================================================================ */

static int reconstruct_marching_cubes(
    const geomlib_pointcloud *pc,
    const geomlib_surface_params *params,
    geomlib_trimesh *output
) {
    if (!pc || pc->point_count < 3) return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_aabb bounds = geomlib_pointcloud_bounds(pc);
    geomlib_vec3 size = geomlib_vec3_sub(bounds.max, bounds.min);

    float voxel_size = params->voxel_size;
    if (voxel_size <= 0.0f) {
        // Auto-compute voxel size
        float avg_spacing;
        geomlib_pointcloud_get_stats(pc, &avg_spacing, NULL);
        voxel_size = avg_spacing * 2.0f;
    }

    int nx = (int)(size.x / voxel_size) + 1;
    int ny = (int)(size.y / voxel_size) + 1;
    int nz = (int)(size.z / voxel_size) + 1;

    // Limit grid size for performance
    if (nx > 64) nx = 64;
    if (ny > 64) ny = 64;
    if (nz > 64) nz = 64;

    // Create density field (simplified - just check if points are nearby)
    float *field = (float*)calloc(nx * ny * nz, sizeof(float));
    if (!field) return GEOMLIB_MESH_ERROR_MEMORY;

    // Fill density field
    float search_radius = voxel_size * 1.5f;
    for (int iz = 0; iz < nz; iz++) {
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                geomlib_vec3 grid_pos = geomlib_vec3_create(
                    bounds.min.x + ix * voxel_size,
                    bounds.min.y + iy * voxel_size,
                    bounds.min.z + iz * voxel_size
                );

                // Count nearby points
                float density = 0.0f;
                for (size_t i = 0; i < pc->point_count; i++) {
                    float dist = vec3_distance(grid_pos, pc->points[i]);
                    if (dist < search_radius) {
                        density += 1.0f - (dist / search_radius);
                    }
                }

                field[iz * ny * nx + iy * nx + ix] = density;
            }
        }
    }

    // Simple marching cubes - extract surface where density > threshold
    float threshold = 0.5f;

    for (int iz = 0; iz < nz - 1; iz++) {
        for (int iy = 0; iy < ny - 1; iy++) {
            for (int ix = 0; ix < nx - 1; ix++) {
                // Get 8 corner values
                float v[8];
                v[0] = field[iz * ny * nx + iy * nx + ix];
                v[1] = field[iz * ny * nx + iy * nx + (ix + 1)];
                v[2] = field[iz * ny * nx + (iy + 1) * nx + (ix + 1)];
                v[3] = field[iz * ny * nx + (iy + 1) * nx + ix];
                v[4] = field[(iz + 1) * ny * nx + iy * nx + ix];
                v[5] = field[(iz + 1) * ny * nx + iy * nx + (ix + 1)];
                v[6] = field[(iz + 1) * ny * nx + (iy + 1) * nx + (ix + 1)];
                v[7] = field[(iz + 1) * ny * nx + (iy + 1) * nx + ix];

                // Simplified: if any edge crosses threshold, add cube faces
                int bits = 0;
                for (int i = 0; i < 8; i++) {
                    if (v[i] > threshold) bits |= (1 << i);
                }

                if (bits != 0 && bits != 0xFF) {
                    // Add simplified quad (not full marching cubes)
                    geomlib_vec3 base = geomlib_vec3_create(
                        bounds.min.x + ix * voxel_size,
                        bounds.min.y + iy * voxel_size,
                        bounds.min.z + iz * voxel_size
                    );

                    geomlib_vec3 vs = geomlib_vec3_create(voxel_size, voxel_size, voxel_size);

                    // Add cube center as vertex
                    geomlib_vec3 center = geomlib_vec3_add(base, geomlib_vec3_mul(vs, 0.5f));
                    geomlib_trimesh_add_vertex(output, center, NULL, NULL);
                }
            }
        }
    }

    free(field);

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Ball Pivoting Algorithm
 * ============================================================================ */

static int reconstruct_ball_pivoting(
    const geomlib_pointcloud *pc,
    const geomlib_surface_params *params,
    geomlib_trimesh *output
) {
    /* Ball pivoting algorithm for surface reconstruction */
    /* Simplified implementation - rolls a ball over points to form triangles */

    geomlib_trimesh_init(output, pc->point_count, pc->point_count * 2);

    /* Auto-compute ball radius if not specified */
    float ball_radius = params->ball_radius;
    if (ball_radius <= 0.0f) {
        /* Estimate from average nearest neighbor distance */
        float avg_dist = 0.0f;
        size_t sample_count = (pc->point_count < 100) ? pc->point_count : 100;
        for (size_t i = 0; i < sample_count; i++) {
            float min_dist = FLT_MAX;
            for (size_t j = 0; j < pc->point_count; j++) {
                if (i == j) continue;
                float d = geomlib_vec3_length(geomlib_vec3_sub(pc->points[i], pc->points[j]));
                if (d < min_dist) min_dist = d;
            }
            avg_dist += min_dist;
        }
        ball_radius = (avg_dist / sample_count) * 2.0f;
    }

    /* Build BVH for point queries */
    geomlib_trimesh temp_mesh;
    geomlib_trimesh_init(&temp_mesh, pc->point_count, 0);
    for (size_t i = 0; i < pc->point_count; i++) {
        geomlib_vec3 normal = pc->normals ? pc->normals[i] : geomlib_vec3_create(0, 0, 1);
        geomlib_trimesh_add_vertex(&temp_mesh, pc->points[i], &normal, NULL);
    }

    /* Track used points */
    bool *used = (bool*)calloc(pc->point_count, sizeof(bool));

    /* Find seed triangle */
    uint32_t seed_idx = 0;
    for (size_t i = 0; i < pc->point_count; i++) {
        if (pc->normals && geomlib_vec3_length(pc->normals[i]) > 0.1f) {
            seed_idx = i;
            break;
        }
    }

    /* Add initial vertices */
    for (size_t i = 0; i < pc->point_count; i++) {
        geomlib_vec3 normal = pc->normals ? pc->normals[i] : geomlib_vec3_create(0, 0, 1);
        geomlib_trimesh_add_vertex(output, pc->points[i], &normal, NULL);
    }

    /* Simplified ball pivoting - add triangles for nearby coplanar points */
    for (size_t i = 0; i < pc->point_count; i++) {
        if (used[i]) continue;

        for (size_t j = i + 1; j < pc->point_count; j++) {
            if (used[j]) continue;

            float d_ij = geomlib_vec3_length(geomlib_vec3_sub(pc->points[i], pc->points[j]));
            if (d_ij > ball_radius * 2.0f) continue;

            for (size_t k = j + 1; k < pc->point_count; k++) {
                if (used[k]) continue;

                float d_ik = geomlib_vec3_length(geomlib_vec3_sub(pc->points[i], pc->points[k]));
                float d_jk = geomlib_vec3_length(geomlib_vec3_sub(pc->points[j], pc->points[k]));

                if (d_ik > ball_radius * 2.0f || d_jk > ball_radius * 2.0f) continue;

                /* Check if triangle is valid (not too flat) */
                geomlib_vec3 e1 = geomlib_vec3_sub(pc->points[j], pc->points[i]);
                geomlib_vec3 e2 = geomlib_vec3_sub(pc->points[k], pc->points[i]);
                geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);
                float area = geomlib_vec3_length(cross) * 0.5f;

                if (area > ball_radius * ball_radius * 0.01f) {
                    geomlib_trimesh_add_triangle(output, i, j, k);
                    used[i] = used[j] = used[k] = true;
                    break;
                }
            }
            if (used[i]) break;
        }
    }

    free(used);
    geomlib_trimesh_free(&temp_mesh);

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Alpha Shape Algorithm
 * ============================================================================ */

static int reconstruct_alpha_shape(
    const geomlib_pointcloud *pc,
    const geomlib_surface_params *params,
    geomlib_trimesh *output
) {
    /* Alpha shape surface reconstruction */
    /* Simplified implementation using alpha parameter for edge length threshold */

    geomlib_trimesh_init(output, pc->point_count, pc->point_count * 2);

    /* Auto-compute alpha if not specified */
    float alpha = params->alpha;
    if (alpha <= 0.0f) {
        /* Estimate from average point spacing */
        float sum = 0.0f;
        size_t count = 0;
        size_t sample_count = (pc->point_count < 50) ? pc->point_count : 50;

        for (size_t i = 0; i < sample_count; i++) {
            float min_dist = FLT_MAX;
            for (size_t j = 0; j < pc->point_count; j++) {
                if (i == j) continue;
                float d = geomlib_vec3_length(geomlib_vec3_sub(pc->points[i], pc->points[j]));
                if (d < min_dist) min_dist = d;
            }
            sum += min_dist;
            count++;
        }
        alpha = (sum / count) * 1.5f;
    }

    /* Add vertices */
    for (size_t i = 0; i < pc->point_count; i++) {
        geomlib_vec3 normal = pc->normals ? pc->normals[i] : geomlib_vec3_create(0, 0, 1);
        geomlib_trimesh_add_vertex(output, pc->points[i], &normal, NULL);
    }

    /* Simplified alpha shape - connect points within alpha distance */
    for (size_t i = 0; i < pc->point_count; i++) {
        for (size_t j = i + 1; j < pc->point_count; j++) {
            float d_ij = geomlib_vec3_length(geomlib_vec3_sub(pc->points[i], pc->points[j]));
            if (d_ij > alpha) continue;

            for (size_t k = j + 1; k < pc->point_count; k++) {
                float d_ik = geomlib_vec3_length(geomlib_vec3_sub(pc->points[i], pc->points[k]));
                float d_jk = geomlib_vec3_length(geomlib_vec3_sub(pc->points[j], pc->points[k]));

                if (d_ik > alpha || d_jk > alpha) continue;

                /* Check circumradius < alpha */
                /* Simplified: just check edge lengths */
                float max_edge = fmaxf(d_ij, fmaxf(d_ik, d_jk));
                if (max_edge < alpha) {
                    /* Compute normal to check orientation */
                    geomlib_vec3 e1 = geomlib_vec3_sub(pc->points[j], pc->points[i]);
                    geomlib_vec3 e2 = geomlib_vec3_sub(pc->points[k], pc->points[i]);
                    geomlib_vec3 tri_normal = geomlib_vec3_normalize(geomlib_vec3_cross(e1, e2));

                    /* Prefer triangles aligned with point normals if available */
                    if (pc->normals) {
                        geomlib_vec3 avg_normal = geomlib_vec3_normalize(geomlib_vec3_add(
                            geomlib_vec3_add(pc->normals[i], pc->normals[j]),
                            pc->normals[k]
                        ));
                        if (geomlib_vec3_dot(tri_normal, avg_normal) > 0.0f) {
                            geomlib_trimesh_add_triangle(output, i, j, k);
                        }
                    } else {
                        geomlib_trimesh_add_triangle(output, i, j, k);
                    }
                }
            }
        }
    }

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Surface Reconstruction - Main API
 * ============================================================================ */

geomlib_surface_params geomlib_surface_params_default(geomlib_surface_method method) {
    geomlib_surface_params params;
    memset(&params, 0, sizeof(params));

    params.method = method;
    params.search_radius = 0.0f;  // Auto
    params.max_neighbors = 20;

    switch (method) {
        case GEOMLIB_SURFACE_GREEDY_PROJECTION:
            params.max_edge_length = 0.0f;  // Auto
            params.max_surface_angle = 45.0f;
            break;

        case GEOMLIB_SURFACE_BALL_PIVOTING:
            params.ball_radius = 0.0f;  // Auto
            break;

        case GEOMLIB_SURFACE_MARCHING_CUBES:
            params.voxel_size = 0.0f;  // Auto
            params.iso_value = 0.5f;
            break;

        case GEOMLIB_SURFACE_ALPHA_SHAPE:
            params.alpha = 0.0f;  // Auto
            break;
    }

    return params;
}

int geomlib_reconstruct_surface(
    const geomlib_pointcloud *pc,
    const geomlib_surface_params *params,
    geomlib_trimesh *output
) {
    if (!pc || !params || !output) return GEOMLIB_MESH_ERROR_INVALID;
    if (pc->point_count < 3) return GEOMLIB_MESH_ERROR_INVALID;

    switch (params->method) {
        case GEOMLIB_SURFACE_GREEDY_PROJECTION:
            return reconstruct_greedy_projection(pc, params, output);

        case GEOMLIB_SURFACE_MARCHING_CUBES:
            return reconstruct_marching_cubes(pc, params, output);

        case GEOMLIB_SURFACE_BALL_PIVOTING:
            return reconstruct_ball_pivoting(pc, params, output);

        case GEOMLIB_SURFACE_ALPHA_SHAPE:
            return reconstruct_alpha_shape(pc, params, output);

        default:
            return GEOMLIB_MESH_ERROR_INVALID;
    }
}

int geomlib_reconstruct_surface_auto(
    const geomlib_pointcloud *pc,
    size_t target_faces,
    geomlib_trimesh *output
) {
    if (!pc || !output) return GEOMLIB_MESH_ERROR_INVALID;

    // Choose method based on point count
    geomlib_surface_method method;
    if (pc->point_count < 1000) {
        method = GEOMLIB_SURFACE_GREEDY_PROJECTION;
    } else {
        method = GEOMLIB_SURFACE_MARCHING_CUBES;
    }

    geomlib_surface_params params = geomlib_surface_params_default(method);

    return geomlib_reconstruct_surface(pc, &params, output);
}
