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
        case GEOMLIB_SURFACE_ALPHA_SHAPE:
            // TODO: implement these methods
            return GEOMLIB_MESH_ERROR_INVALID;

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
