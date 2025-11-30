/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Monte Carlo Geometry Processing Implementation
 */

#include "geomlib/geomlib_montecarlo.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Triangle Surface Sampling
 * ============================================================================ */

geomlib_vec3 geomlib_sample_triangle_uniform(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_pcg32 *rng
) {
    /* Standard uniform triangle sampling using barycentric coordinates
     * Reference: pbr-book.org, Matt Pharr's blog
     * Formula: b0 = 1 - sqrt(u), b1 = v * sqrt(u), b2 = 1 - b0 - b1
     */
    float u = geomlib_pcg32_nextf(rng);
    float v = geomlib_pcg32_nextf(rng);

    float sqrt_u = sqrtf(u);
    float b0 = 1.0f - sqrt_u;
    float b1 = v * sqrt_u;
    float b2 = 1.0f - b0 - b1;

    /* Compute point using barycentric coordinates */
    geomlib_vec3 result;
    result.x = b0 * v0.x + b1 * v1.x + b2 * v2.x;
    result.y = b0 * v0.y + b1 * v1.y + b2 * v2.y;
    result.z = b0 * v0.z + b1 * v1.z + b2 * v2.z;

    return result;
}

geomlib_vec3 geomlib_sample_triangle_barycentric(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    geomlib_pcg32 *rng,
    geomlib_vec3 *bary_out
) {
    float u = geomlib_pcg32_nextf(rng);
    float v = geomlib_pcg32_nextf(rng);

    float sqrt_u = sqrtf(u);
    float b0 = 1.0f - sqrt_u;
    float b1 = v * sqrt_u;
    float b2 = 1.0f - b0 - b1;

    if (bary_out) {
        bary_out->x = b0;
        bary_out->y = b1;
        bary_out->z = b2;
    }

    geomlib_vec3 result;
    result.x = b0 * v0.x + b1 * v1.x + b2 * v2.x;
    result.y = b0 * v0.y + b1 * v1.y + b2 * v2.y;
    result.z = b0 * v0.z + b1 * v1.z + b2 * v2.z;

    return result;
}

float geomlib_triangle_sample_pdf(
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2
) {
    /* PDF for uniform sampling = 1 / area */
    geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);
    float area = geomlib_vec3_length(cross) * 0.5f;

    return (area > 0.0f) ? (1.0f / area) : 0.0f;
}

/* ============================================================================
 * Mesh Surface Sampling
 * ============================================================================ */

int geomlib_mesh_build_sample_distribution(
    const geomlib_trimesh *mesh,
    geomlib_mesh_sample_distribution *dist
) {
    if (!mesh || !dist || mesh->triangle_count == 0) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    dist->triangle_count = mesh->triangle_count;
    dist->triangle_areas = (float *)malloc(sizeof(float) * mesh->triangle_count);
    dist->cdf = (float *)malloc(sizeof(float) * mesh->triangle_count);

    if (!dist->triangle_areas || !dist->cdf) {
        free(dist->triangle_areas);
        free(dist->cdf);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    /* Compute triangle areas and build CDF */
    dist->total_area = 0.0f;

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t i0 = mesh->indices[i * 3 + 0];
        uint32_t i1 = mesh->indices[i * 3 + 1];
        uint32_t i2 = mesh->indices[i * 3 + 2];

        geomlib_vec3 v0 = mesh->vertices[i0];
        geomlib_vec3 v1 = mesh->vertices[i1];
        geomlib_vec3 v2 = mesh->vertices[i2];

        geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
        geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
        geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);
        float area = geomlib_vec3_length(cross) * 0.5f;

        dist->triangle_areas[i] = area;
        dist->total_area += area;
        dist->cdf[i] = dist->total_area;
    }

    /* Normalize CDF to [0, 1] */
    if (dist->total_area > 0.0f) {
        for (size_t i = 0; i < mesh->triangle_count; i++) {
            dist->cdf[i] /= dist->total_area;
        }
    }

    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_mesh_free_sample_distribution(
    geomlib_mesh_sample_distribution *dist
) {
    if (dist) {
        free(dist->triangle_areas);
        free(dist->cdf);
        dist->triangle_areas = NULL;
        dist->cdf = NULL;
        dist->triangle_count = 0;
        dist->total_area = 0.0f;
    }
}

static size_t binary_search_cdf(const float *cdf, size_t count, float value) {
    /* Binary search to find index where cdf[i-1] <= value < cdf[i] */
    size_t left = 0;
    size_t right = count;

    while (left < right) {
        size_t mid = left + (right - left) / 2;
        if (cdf[mid] <= value) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    return (left < count) ? left : count - 1;
}

geomlib_vec3 geomlib_sample_mesh_surface(
    const geomlib_trimesh *mesh,
    const geomlib_mesh_sample_distribution *dist,
    geomlib_pcg32 *rng,
    size_t *triangle_idx_out
) {
    /* Sample triangle using CDF (binary search) */
    float u = geomlib_pcg32_nextf(rng);
    size_t tri_idx = binary_search_cdf(dist->cdf, dist->triangle_count, u);

    if (triangle_idx_out) {
        *triangle_idx_out = tri_idx;
    }

    /* Get triangle vertices */
    uint32_t i0 = mesh->indices[tri_idx * 3 + 0];
    uint32_t i1 = mesh->indices[tri_idx * 3 + 1];
    uint32_t i2 = mesh->indices[tri_idx * 3 + 2];

    geomlib_vec3 v0 = mesh->vertices[i0];
    geomlib_vec3 v1 = mesh->vertices[i1];
    geomlib_vec3 v2 = mesh->vertices[i2];

    /* Sample point on triangle */
    return geomlib_sample_triangle_uniform(v0, v1, v2, rng);
}

geomlib_vec3 geomlib_sample_mesh_surface_with_normal(
    const geomlib_trimesh *mesh,
    const geomlib_mesh_sample_distribution *dist,
    geomlib_pcg32 *rng,
    geomlib_vec3 *normal_out,
    size_t *triangle_idx_out
) {
    geomlib_vec3 bary;

    /* Sample triangle using CDF */
    float u = geomlib_pcg32_nextf(rng);
    size_t tri_idx = binary_search_cdf(dist->cdf, dist->triangle_count, u);

    if (triangle_idx_out) {
        *triangle_idx_out = tri_idx;
    }

    /* Get triangle vertices */
    uint32_t i0 = mesh->indices[tri_idx * 3 + 0];
    uint32_t i1 = mesh->indices[tri_idx * 3 + 1];
    uint32_t i2 = mesh->indices[tri_idx * 3 + 2];

    geomlib_vec3 v0 = mesh->vertices[i0];
    geomlib_vec3 v1 = mesh->vertices[i1];
    geomlib_vec3 v2 = mesh->vertices[i2];

    /* Sample point on triangle with barycentric coordinates */
    geomlib_vec3 point = geomlib_sample_triangle_barycentric(v0, v1, v2, rng, &bary);

    /* Interpolate normal if available */
    if (normal_out) {
        if (mesh->normals) {
            geomlib_vec3 n0 = mesh->normals[i0];
            geomlib_vec3 n1 = mesh->normals[i1];
            geomlib_vec3 n2 = mesh->normals[i2];

            normal_out->x = bary.x * n0.x + bary.y * n1.x + bary.z * n2.x;
            normal_out->y = bary.x * n0.y + bary.y * n1.y + bary.z * n2.y;
            normal_out->z = bary.x * n0.z + bary.y * n1.z + bary.z * n2.z;

            *normal_out = geomlib_vec3_normalize(*normal_out);
        } else {
            /* Compute face normal */
            geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
            geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
            *normal_out = geomlib_vec3_normalize(geomlib_vec3_cross(e1, e2));
        }
    }

    return point;
}

/* ============================================================================
 * Monte Carlo Integration Helpers
 * ============================================================================ */

double geomlib_mc_estimate_surface_integral(
    const geomlib_vec3 *samples,
    const float *values,
    const float *pdf,
    size_t num_samples,
    float area
) {
    if (num_samples == 0) {
        return 0.0;
    }

    double sum = 0.0;

    if (pdf) {
        /* General case: (1/N) * Σ f(x_i) / pdf(x_i) */
        for (size_t i = 0; i < num_samples; i++) {
            if (pdf[i] > 0.0f) {
                sum += values[i] / pdf[i];
            }
        }
        return sum / num_samples;
    } else {
        /* Uniform sampling: (area/N) * Σ f(x_i) */
        for (size_t i = 0; i < num_samples; i++) {
            sum += values[i];
        }
        return (sum / num_samples) * area;
    }
}

double geomlib_mc_estimate_variance(
    const float *values,
    const float *pdf,
    size_t num_samples,
    double mean
) {
    if (num_samples <= 1) {
        return 0.0;
    }

    double sum_sq_diff = 0.0;

    if (pdf) {
        for (size_t i = 0; i < num_samples; i++) {
            if (pdf[i] > 0.0f) {
                double sample = values[i] / pdf[i];
                double diff = sample - mean;
                sum_sq_diff += diff * diff;
            }
        }
    } else {
        for (size_t i = 0; i < num_samples; i++) {
            double diff = values[i] - mean;
            sum_sq_diff += diff * diff;
        }
    }

    return sum_sq_diff / (num_samples - 1);
}

/* ============================================================================
 * Closest Point Queries
 * ============================================================================ */

geomlib_vec3 geomlib_closest_point_triangle(
    geomlib_vec3 query,
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2
) {
    /* Based on Real-Time Collision Detection by Christer Ericson */
    geomlib_vec3 edge0 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 edge1 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 v0_to_p = geomlib_vec3_sub(query, v0);

    float a = geomlib_vec3_dot(edge0, edge0);
    float b = geomlib_vec3_dot(edge0, edge1);
    float c = geomlib_vec3_dot(edge1, edge1);
    float d = geomlib_vec3_dot(edge0, v0_to_p);
    float e = geomlib_vec3_dot(edge1, v0_to_p);

    float det = a * c - b * b;
    float s = c * d - b * e;
    float t = a * e - b * d;

    if (s + t <= det) {
        if (s < 0.0f) {
            if (t < 0.0f) {
                /* Region 4: closest to v0 */
                s = 0.0f;
                t = 0.0f;
            } else {
                /* Region 3: closest to edge v0-v2 */
                s = 0.0f;
                t = (e >= 0.0f) ? 0.0f : ((e <= -c) ? 1.0f : -e / c);
            }
        } else if (t < 0.0f) {
            /* Region 5: closest to edge v0-v1 */
            t = 0.0f;
            s = (d >= 0.0f) ? 0.0f : ((d <= -a) ? 1.0f : -d / a);
        } else {
            /* Region 0: interior point */
            float inv_det = 1.0f / det;
            s *= inv_det;
            t *= inv_det;
        }
    } else {
        if (s < 0.0f) {
            /* Region 2: closest to edge v1-v2 */
            float tmp0 = b + d;
            float tmp1 = c + e;
            if (tmp1 > tmp0) {
                float numer = tmp1 - tmp0;
                float denom = a - 2.0f * b + c;
                s = (numer >= denom) ? 1.0f : numer / denom;
                t = 1.0f - s;
            } else {
                s = 0.0f;
                t = (tmp1 <= 0.0f) ? 1.0f : ((e >= 0.0f) ? 0.0f : -e / c);
            }
        } else if (t < 0.0f) {
            /* Region 6: closest to edge v0-v1 */
            t = 0.0f;
            s = (d >= 0.0f) ? 0.0f : ((d <= -a) ? 1.0f : -d / a);
        } else {
            /* Region 1: closest to edge v1-v2 */
            float numer = (c + e) - (b + d);
            if (numer <= 0.0f) {
                s = 0.0f;
            } else {
                float denom = a - 2.0f * b + c;
                s = (numer >= denom) ? 1.0f : numer / denom;
            }
            t = 1.0f - s;
        }
    }

    /* Compute closest point */
    geomlib_vec3 result = v0;
    result = geomlib_vec3_add(result, geomlib_vec3_mul(edge0, s));
    result = geomlib_vec3_add(result, geomlib_vec3_mul(edge1, t));

    return result;
}

geomlib_vec3 geomlib_closest_point_mesh(
    geomlib_vec3 query,
    const geomlib_trimesh *mesh,
    float *distance_out
) {
    geomlib_vec3 closest_point = {0, 0, 0};
    float min_dist_sq = FLT_MAX;

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t i0 = mesh->indices[i * 3 + 0];
        uint32_t i1 = mesh->indices[i * 3 + 1];
        uint32_t i2 = mesh->indices[i * 3 + 2];

        geomlib_vec3 v0 = mesh->vertices[i0];
        geomlib_vec3 v1 = mesh->vertices[i1];
        geomlib_vec3 v2 = mesh->vertices[i2];

        geomlib_vec3 pt = geomlib_closest_point_triangle(query, v0, v1, v2);
        geomlib_vec3 diff = geomlib_vec3_sub(pt, query);
        float dist_sq = geomlib_vec3_dot(diff, diff);

        if (dist_sq < min_dist_sq) {
            min_dist_sq = dist_sq;
            closest_point = pt;
        }
    }

    if (distance_out) {
        *distance_out = sqrtf(min_dist_sq);
    }

    return closest_point;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

geomlib_vec3 geomlib_sample_sphere_uniform(geomlib_pcg32 *rng) {
    /* Uniform sampling on unit sphere using spherical coordinates */
    float u = geomlib_pcg32_nextf(rng);
    float v = geomlib_pcg32_nextf(rng);

    float theta = 2.0f * M_PI * u;
    float phi = acosf(2.0f * v - 1.0f);

    float sin_phi = sinf(phi);

    geomlib_vec3 result;
    result.x = sin_phi * cosf(theta);
    result.y = sin_phi * sinf(theta);
    result.z = cosf(phi);

    return result;
}

geomlib_vec3 geomlib_sample_sphere_volume(geomlib_pcg32 *rng) {
    /* Rejection sampling in unit cube */
    geomlib_vec3 p;
    float dist_sq;

    do {
        p.x = geomlib_pcg32_range_float(rng, -1.0f, 1.0f);
        p.y = geomlib_pcg32_range_float(rng, -1.0f, 1.0f);
        p.z = geomlib_pcg32_range_float(rng, -1.0f, 1.0f);
        dist_sq = p.x * p.x + p.y * p.y + p.z * p.z;
    } while (dist_sq > 1.0f);

    return p;
}

geomlib_vec3 geomlib_sample_hemisphere_uniform(geomlib_vec3 normal, geomlib_pcg32 *rng) {
    /* Sample uniform direction in hemisphere */
    float u = geomlib_pcg32_nextf(rng);
    float v = geomlib_pcg32_nextf(rng);

    float theta = 2.0f * M_PI * u;
    float phi = acosf(v);

    float sin_phi = sinf(phi);

    /* Local direction */
    geomlib_vec3 local_dir;
    local_dir.x = sin_phi * cosf(theta);
    local_dir.y = sin_phi * sinf(theta);
    local_dir.z = cosf(phi);

    /* Build orthonormal basis from normal */
    geomlib_vec3 tangent, bitangent;

    if (fabsf(normal.x) > fabsf(normal.y)) {
        float inv_len = 1.0f / sqrtf(normal.x * normal.x + normal.z * normal.z);
        tangent.x = -normal.z * inv_len;
        tangent.y = 0.0f;
        tangent.z = normal.x * inv_len;
    } else {
        float inv_len = 1.0f / sqrtf(normal.y * normal.y + normal.z * normal.z);
        tangent.x = 0.0f;
        tangent.y = normal.z * inv_len;
        tangent.z = -normal.y * inv_len;
    }

    bitangent = geomlib_vec3_cross(normal, tangent);

    /* Transform to world space */
    geomlib_vec3 result;
    result.x = local_dir.x * tangent.x + local_dir.y * bitangent.x + local_dir.z * normal.x;
    result.y = local_dir.x * tangent.y + local_dir.y * bitangent.y + local_dir.z * normal.y;
    result.z = local_dir.x * tangent.z + local_dir.y * bitangent.z + local_dir.z * normal.z;

    return result;
}

geomlib_vec3 geomlib_sample_hemisphere_cosine(geomlib_vec3 normal, geomlib_pcg32 *rng) {
    /* Cosine-weighted hemisphere sampling (Malley's method) */
    float u = geomlib_pcg32_nextf(rng);
    float v = geomlib_pcg32_nextf(rng);

    float r = sqrtf(u);
    float theta = 2.0f * M_PI * v;

    /* Local direction on unit disk, then project to hemisphere */
    geomlib_vec3 local_dir;
    local_dir.x = r * cosf(theta);
    local_dir.y = r * sinf(theta);
    local_dir.z = sqrtf(fmaxf(0.0f, 1.0f - u));

    /* Build orthonormal basis */
    geomlib_vec3 tangent, bitangent;

    if (fabsf(normal.x) > fabsf(normal.y)) {
        float inv_len = 1.0f / sqrtf(normal.x * normal.x + normal.z * normal.z);
        tangent.x = -normal.z * inv_len;
        tangent.y = 0.0f;
        tangent.z = normal.x * inv_len;
    } else {
        float inv_len = 1.0f / sqrtf(normal.y * normal.y + normal.z * normal.z);
        tangent.x = 0.0f;
        tangent.y = normal.z * inv_len;
        tangent.z = -normal.y * inv_len;
    }

    bitangent = geomlib_vec3_cross(normal, tangent);

    /* Transform to world space */
    geomlib_vec3 result;
    result.x = local_dir.x * tangent.x + local_dir.y * bitangent.x + local_dir.z * normal.x;
    result.y = local_dir.x * tangent.y + local_dir.y * bitangent.y + local_dir.z * normal.y;
    result.z = local_dir.x * tangent.z + local_dir.y * bitangent.z + local_dir.z * normal.z;

    return result;
}
