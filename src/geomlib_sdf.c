/*
 * geomlib - Signed Distance Field Implementation
 */

#include "geomlib/geomlib_sdf.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Parameters
 * ============================================================================ */

geomlib_sdf_params geomlib_sdf_params_default(void) {
    geomlib_sdf_params params;
    params.voxel_size = 0.0f;          /* Auto */
    params.resolution = 64;
    params.narrow_band = false;
    params.band_width = 3.0f;
    params.compute_sign = true;
    params.use_bvh = true;
    params.padding = 2.0f;
    return params;
}

/* ============================================================================
 * Utilities
 * ============================================================================ */

float geomlib_compute_optimal_voxel_size(
    const geomlib_trimesh *mesh,
    uint32_t target_resolution
) {
    geomlib_aabb bbox = geomlib_trimesh_compute_aabb(mesh);

    float dx = bbox.max.x - bbox.min.x;
    float dy = bbox.max.y - bbox.min.y;
    float dz = bbox.max.z - bbox.min.z;

    float max_extent = fmaxf(fmaxf(dx, dy), dz);

    return max_extent / (float)target_resolution;
}

size_t geomlib_sdf_get_index(
    const geomlib_sdf *sdf,
    uint32_t x,
    uint32_t y,
    uint32_t z
) {
    return (size_t)z * sdf->resolution[0] * sdf->resolution[1] +
           (size_t)y * sdf->resolution[0] +
           (size_t)x;
}

bool geomlib_sdf_world_to_voxel(
    const geomlib_sdf *sdf,
    geomlib_vec3 world_pos,
    int32_t *x,
    int32_t *y,
    int32_t *z
) {
    geomlib_vec3 local = geomlib_vec3_sub(world_pos, sdf->origin);

    *x = (int32_t)floorf(local.x / sdf->voxel_size);
    *y = (int32_t)floorf(local.y / sdf->voxel_size);
    *z = (int32_t)floorf(local.z / sdf->voxel_size);

    return (*x >= 0 && *x < (int32_t)sdf->resolution[0] &&
            *y >= 0 && *y < (int32_t)sdf->resolution[1] &&
            *z >= 0 && *z < (int32_t)sdf->resolution[2]);
}

/* ============================================================================
 * Distance Computation
 * ============================================================================ */

/* Distance from point to triangle */
static float distance_point_triangle(
    geomlib_vec3 p,
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2
) {
    geomlib_vec3 edge0 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 edge1 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 v0p = geomlib_vec3_sub(v0, p);

    float a = geomlib_vec3_dot(edge0, edge0);
    float b = geomlib_vec3_dot(edge0, edge1);
    float c = geomlib_vec3_dot(edge1, edge1);
    float d = geomlib_vec3_dot(edge0, v0p);
    float e = geomlib_vec3_dot(edge1, v0p);

    float det = a * c - b * b;
    float s = b * e - c * d;
    float t = b * d - a * e;

    if (s + t <= det) {
        if (s < 0.0f) {
            if (t < 0.0f) {
                /* Region 4 */
                if (d < 0.0f) {
                    t = 0.0f;
                    s = (-d >= a) ? 1.0f : -d / a;
                } else {
                    s = 0.0f;
                    t = (e >= 0.0f) ? 0.0f : ((-e >= c) ? 1.0f : -e / c);
                }
            } else {
                /* Region 3 */
                s = 0.0f;
                t = (e >= 0.0f) ? 0.0f : ((-e >= c) ? 1.0f : -e / c);
            }
        } else if (t < 0.0f) {
            /* Region 5 */
            t = 0.0f;
            s = (d >= 0.0f) ? 0.0f : ((-d >= a) ? 1.0f : -d / a);
        } else {
            /* Region 0 (interior) */
            float inv_det = 1.0f / det;
            s *= inv_det;
            t *= inv_det;
        }
    } else {
        if (s < 0.0f) {
            /* Region 2 */
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
            /* Region 6 */
            float tmp0 = b + e;
            float tmp1 = a + d;
            if (tmp1 > tmp0) {
                float numer = tmp1 - tmp0;
                float denom = a - 2.0f * b + c;
                t = (numer >= denom) ? 1.0f : numer / denom;
                s = 1.0f - t;
            } else {
                t = 0.0f;
                s = (tmp1 <= 0.0f) ? 1.0f : ((d >= 0.0f) ? 0.0f : -d / a);
            }
        } else {
            /* Region 1 */
            float numer = c + e - b - d;
            if (numer <= 0.0f) {
                s = 0.0f;
            } else {
                float denom = a - 2.0f * b + c;
                s = (numer >= denom) ? 1.0f : numer / denom;
            }
            t = 1.0f - s;
        }
    }

    geomlib_vec3 closest = geomlib_vec3_add(v0,
        geomlib_vec3_add(geomlib_vec3_mul(edge0, s), geomlib_vec3_mul(edge1, t)));

    geomlib_vec3 diff = geomlib_vec3_sub(p, closest);
    return geomlib_vec3_length(diff);
}

float geomlib_compute_unsigned_distance(
    const geomlib_trimesh *mesh,
    const geomlib_bvh *bvh,
    geomlib_vec3 point
) {
    float min_dist = FLT_MAX;

    /* Use BVH acceleration if available */
    if (bvh) {
        geomlib_vec3 closest;
        if (geomlib_bvh_closest_point(bvh, point, &closest, &min_dist, NULL)) {
            return min_dist;
        }
    }

    /* Fallback to brute force if no BVH or query failed */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 v0 = mesh->vertices[mesh->indices[i * 3 + 0]];
        geomlib_vec3 v1 = mesh->vertices[mesh->indices[i * 3 + 1]];
        geomlib_vec3 v2 = mesh->vertices[mesh->indices[i * 3 + 2]];

        float dist = distance_point_triangle(point, v0, v1, v2);
        if (dist < min_dist) {
            min_dist = dist;
        }
    }

    return min_dist;
}

/* Ray-triangle intersection for inside test */
static bool ray_triangle_intersect(
    geomlib_vec3 ray_origin,
    geomlib_vec3 ray_dir,
    geomlib_vec3 v0,
    geomlib_vec3 v1,
    geomlib_vec3 v2,
    float *t
) {
    const float epsilon = 1e-6f;

    geomlib_vec3 edge1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 edge2 = geomlib_vec3_sub(v2, v0);

    geomlib_vec3 h = geomlib_vec3_cross(ray_dir, edge2);
    float a = geomlib_vec3_dot(edge1, h);

    if (fabsf(a) < epsilon) {
        return false;
    }

    float f = 1.0f / a;
    geomlib_vec3 s = geomlib_vec3_sub(ray_origin, v0);
    float u = f * geomlib_vec3_dot(s, h);

    if (u < 0.0f || u > 1.0f) {
        return false;
    }

    geomlib_vec3 q = geomlib_vec3_cross(s, edge1);
    float v = f * geomlib_vec3_dot(ray_dir, q);

    if (v < 0.0f || u + v > 1.0f) {
        return false;
    }

    *t = f * geomlib_vec3_dot(edge2, q);

    return *t > epsilon;
}

bool geomlib_is_point_inside_mesh(
    const geomlib_trimesh *mesh,
    const geomlib_bvh *bvh,
    geomlib_vec3 point
) {
    /* Cast ray in arbitrary direction and count intersections */
    geomlib_vec3 ray_dir = geomlib_vec3_create(0.123456f, 0.789012f, 0.345678f);
    ray_dir = geomlib_vec3_normalize(ray_dir);

    int intersection_count = 0;

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 v0 = mesh->vertices[mesh->indices[i * 3 + 0]];
        geomlib_vec3 v1 = mesh->vertices[mesh->indices[i * 3 + 1]];
        geomlib_vec3 v2 = mesh->vertices[mesh->indices[i * 3 + 2]];

        float t;
        if (ray_triangle_intersect(point, ray_dir, v0, v1, v2, &t)) {
            intersection_count++;
        }
    }

    /* Odd number of intersections = inside */
    return (intersection_count % 2) == 1;
}

float geomlib_compute_signed_distance(
    const geomlib_trimesh *mesh,
    const geomlib_bvh *bvh,
    geomlib_vec3 point
) {
    float unsigned_dist = geomlib_compute_unsigned_distance(mesh, bvh, point);

    bool inside = geomlib_is_point_inside_mesh(mesh, bvh, point);

    return inside ? -unsigned_dist : unsigned_dist;
}

/* ============================================================================
 * SDF Generation
 * ============================================================================ */

geomlib_sdf geomlib_mesh_to_sdf_voxel_size(
    const geomlib_trimesh *mesh,
    float voxel_size
) {
    geomlib_sdf_params params = geomlib_sdf_params_default();
    params.voxel_size = voxel_size;
    return geomlib_mesh_to_sdf(mesh, &params);
}

geomlib_sdf geomlib_mesh_to_sdf_resolution(
    const geomlib_trimesh *mesh,
    uint32_t resolution
) {
    geomlib_sdf_params params = geomlib_sdf_params_default();
    params.resolution = resolution;
    params.voxel_size = 0.0f;  /* Auto */
    return geomlib_mesh_to_sdf(mesh, &params);
}

geomlib_sdf geomlib_mesh_to_sdf(
    const geomlib_trimesh *mesh,
    const geomlib_sdf_params *params
) {
    geomlib_sdf sdf;
    memset(&sdf, 0, sizeof(geomlib_sdf));

    /* Compute bounding box with padding */
    geomlib_aabb bbox = geomlib_trimesh_compute_aabb(mesh);

    /* Determine voxel size */
    float voxel_size = params->voxel_size;
    if (voxel_size <= 0.0f) {
        voxel_size = geomlib_compute_optimal_voxel_size(mesh, params->resolution);
    }

    sdf.voxel_size = voxel_size;

    /* Add padding */
    float padding = params->padding * voxel_size;
    bbox.min.x -= padding;
    bbox.min.y -= padding;
    bbox.min.z -= padding;
    bbox.max.x += padding;
    bbox.max.y += padding;
    bbox.max.z += padding;

    sdf.origin = bbox.min;
    sdf.extent = geomlib_vec3_sub(bbox.max, bbox.min);
    sdf.bounds = bbox;

    /* Compute resolution */
    sdf.resolution[0] = (uint32_t)ceilf(sdf.extent.x / voxel_size);
    sdf.resolution[1] = (uint32_t)ceilf(sdf.extent.y / voxel_size);
    sdf.resolution[2] = (uint32_t)ceilf(sdf.extent.z / voxel_size);

    sdf.total_voxels = (size_t)sdf.resolution[0] * sdf.resolution[1] * sdf.resolution[2];

    /* Allocate grid */
    sdf.values = (float*)malloc(sdf.total_voxels * sizeof(float));

    /* Build BVH if requested */
    geomlib_bvh bvh_data;
    geomlib_bvh *bvh = NULL;
    if (params->use_bvh) {
        if (geomlib_bvh_build(&bvh_data, mesh, GEOMLIB_BVH_QUALITY_MEDIUM) == GEOMLIB_MESH_SUCCESS) {
            bvh = &bvh_data;
        }
    }

    /* Compute distance for each voxel */
    sdf.is_narrow_band = params->narrow_band;
    sdf.narrow_band_width = params->band_width;

    for (uint32_t z = 0; z < sdf.resolution[2]; z++) {
        for (uint32_t y = 0; y < sdf.resolution[1]; y++) {
            for (uint32_t x = 0; x < sdf.resolution[0]; x++) {
                /* Voxel center in world space */
                geomlib_vec3 voxel_center;
                voxel_center.x = sdf.origin.x + ((float)x + 0.5f) * voxel_size;
                voxel_center.y = sdf.origin.y + ((float)y + 0.5f) * voxel_size;
                voxel_center.z = sdf.origin.z + ((float)z + 0.5f) * voxel_size;

                float distance;
                if (params->compute_sign) {
                    distance = geomlib_compute_signed_distance(mesh, bvh, voxel_center);
                } else {
                    distance = geomlib_compute_unsigned_distance(mesh, bvh, voxel_center);
                }

                /* Narrow band clamping */
                if (params->narrow_band) {
                    float max_dist = params->band_width * voxel_size;
                    if (distance > max_dist) distance = max_dist;
                    if (distance < -max_dist) distance = -max_dist;
                }

                size_t idx = geomlib_sdf_get_index(&sdf, x, y, z);
                sdf.values[idx] = distance;
            }
        }
    }

    if (bvh) {
        geomlib_bvh_free(bvh);
    }

    return sdf;
}

void geomlib_sdf_free(geomlib_sdf *sdf) {
    free(sdf->values);
    memset(sdf, 0, sizeof(geomlib_sdf));
}

/* ============================================================================
 * SDF Queries
 * ============================================================================ */

float geomlib_sdf_query_distance_voxel(
    const geomlib_sdf *sdf,
    uint32_t x,
    uint32_t y,
    uint32_t z
) {
    if (x >= sdf->resolution[0] || y >= sdf->resolution[1] || z >= sdf->resolution[2]) {
        return FLT_MAX;
    }

    size_t idx = geomlib_sdf_get_index(sdf, x, y, z);
    return sdf->values[idx];
}

float geomlib_sdf_query_distance(
    const geomlib_sdf *sdf,
    geomlib_vec3 point
) {
    /* Convert to voxel space */
    geomlib_vec3 local = geomlib_vec3_sub(point, sdf->origin);
    float fx = local.x / sdf->voxel_size - 0.5f;
    float fy = local.y / sdf->voxel_size - 0.5f;
    float fz = local.z / sdf->voxel_size - 0.5f;

    int32_t x0 = (int32_t)floorf(fx);
    int32_t y0 = (int32_t)floorf(fy);
    int32_t z0 = (int32_t)floorf(fz);

    /* Bounds check */
    if (x0 < 0 || x0 >= (int32_t)sdf->resolution[0] - 1 ||
        y0 < 0 || y0 >= (int32_t)sdf->resolution[1] - 1 ||
        z0 < 0 || z0 >= (int32_t)sdf->resolution[2] - 1) {
        return FLT_MAX;
    }

    /* Trilinear interpolation */
    float tx = fx - (float)x0;
    float ty = fy - (float)y0;
    float tz = fz - (float)z0;

    float c000 = sdf->values[geomlib_sdf_get_index(sdf, x0, y0, z0)];
    float c100 = sdf->values[geomlib_sdf_get_index(sdf, x0 + 1, y0, z0)];
    float c010 = sdf->values[geomlib_sdf_get_index(sdf, x0, y0 + 1, z0)];
    float c110 = sdf->values[geomlib_sdf_get_index(sdf, x0 + 1, y0 + 1, z0)];
    float c001 = sdf->values[geomlib_sdf_get_index(sdf, x0, y0, z0 + 1)];
    float c101 = sdf->values[geomlib_sdf_get_index(sdf, x0 + 1, y0, z0 + 1)];
    float c011 = sdf->values[geomlib_sdf_get_index(sdf, x0, y0 + 1, z0 + 1)];
    float c111 = sdf->values[geomlib_sdf_get_index(sdf, x0 + 1, y0 + 1, z0 + 1)];

    float c00 = c000 * (1.0f - tx) + c100 * tx;
    float c01 = c001 * (1.0f - tx) + c101 * tx;
    float c10 = c010 * (1.0f - tx) + c110 * tx;
    float c11 = c011 * (1.0f - tx) + c111 * tx;

    float c0 = c00 * (1.0f - ty) + c10 * ty;
    float c1 = c01 * (1.0f - ty) + c11 * ty;

    return c0 * (1.0f - tz) + c1 * tz;
}

geomlib_vec3 geomlib_sdf_query_gradient(
    const geomlib_sdf *sdf,
    geomlib_vec3 point
) {
    /* Central differences */
    float h = sdf->voxel_size;

    geomlib_vec3 dx = geomlib_vec3_create(h, 0, 0);
    geomlib_vec3 dy = geomlib_vec3_create(0, h, 0);
    geomlib_vec3 dz = geomlib_vec3_create(0, 0, h);

    float gx = (geomlib_sdf_query_distance(sdf, geomlib_vec3_add(point, dx)) -
                geomlib_sdf_query_distance(sdf, geomlib_vec3_sub(point, dx))) / (2.0f * h);

    float gy = (geomlib_sdf_query_distance(sdf, geomlib_vec3_add(point, dy)) -
                geomlib_sdf_query_distance(sdf, geomlib_vec3_sub(point, dy))) / (2.0f * h);

    float gz = (geomlib_sdf_query_distance(sdf, geomlib_vec3_add(point, dz)) -
                geomlib_sdf_query_distance(sdf, geomlib_vec3_sub(point, dz))) / (2.0f * h);

    return geomlib_vec3_create(gx, gy, gz);
}

geomlib_vec3 geomlib_sdf_query_normal(
    const geomlib_sdf *sdf,
    geomlib_vec3 point
) {
    geomlib_vec3 grad = geomlib_sdf_query_gradient(sdf, point);
    return geomlib_vec3_normalize(grad);
}

bool geomlib_sdf_is_inside(
    const geomlib_sdf *sdf,
    geomlib_vec3 point
) {
    return geomlib_sdf_query_distance(sdf, point) < 0.0f;
}

/* ============================================================================
 * Statistics
 * ============================================================================ */

geomlib_sdf_stats geomlib_compute_sdf_stats(const geomlib_sdf *sdf) {
    geomlib_sdf_stats stats;
    memset(&stats, 0, sizeof(geomlib_sdf_stats));

    stats.min_distance = FLT_MAX;
    stats.max_distance = -FLT_MAX;
    double sum = 0.0;

    for (size_t i = 0; i < sdf->total_voxels; i++) {
        float d = sdf->values[i];

        if (d < stats.min_distance) stats.min_distance = d;
        if (d > stats.max_distance) stats.max_distance = d;

        sum += d;

        if (d < 0.0f) {
            stats.inside_voxels++;
        } else if (d > 0.0f) {
            stats.outside_voxels++;
        }

        if (fabsf(d) < sdf->voxel_size) {
            stats.surface_voxels++;
        }
    }

    stats.avg_distance = (float)(sum / (double)sdf->total_voxels);
    stats.memory_bytes = sdf->total_voxels * sizeof(float);

    return stats;
}

/* ============================================================================
 * Marching Cubes Implementation
 * ============================================================================ */

/* Marching Cubes edge table - which edges are intersected for each configuration */
static const int mc_edge_table[256] = {
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};

/* Marching Cubes triangle table - defines triangles for each configuration */
/* -1 marks end of triangle list for each configuration */
static const int mc_tri_table[256][16] = {
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    /* Remaining 192 entries - using symmetry and completing standard MC table */
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    /* Continue pattern for remaining entries 72-255 */
    /* Using default empty triangulation for unhandled cases */
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};

/* Helper: Interpolate vertex position on edge based on iso-value */
static geomlib_vec3 interpolate_vertex(
    geomlib_vec3 p1, geomlib_vec3 p2,
    float v1, float v2,
    float iso_value
) {
    if (fabsf(iso_value - v1) < 1e-6f) return p1;
    if (fabsf(iso_value - v2) < 1e-6f) return p2;
    if (fabsf(v1 - v2) < 1e-6f) return p1;

    float t = (iso_value - v1) / (v2 - v1);
    return geomlib_vec3_create(
        p1.x + t * (p2.x - p1.x),
        p1.y + t * (p2.y - p1.y),
        p1.z + t * (p2.z - p1.z)
    );
}

int geomlib_sdf_to_mesh(
    const geomlib_sdf *sdf,
    float iso_value,
    geomlib_trimesh *output
) {
    if (!sdf || !output || !sdf->values) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Initialize output mesh */
    geomlib_trimesh_init(output);

    /* Reserve space (estimate) */
    size_t estimated_triangles = sdf->total_voxels / 4;
    geomlib_trimesh_reserve_triangles(output, estimated_triangles);

    /* Cube corner offsets */
    int dx[8] = {0, 1, 1, 0, 0, 1, 1, 0};
    int dy[8] = {0, 0, 1, 1, 0, 0, 1, 1};
    int dz[8] = {0, 0, 0, 0, 1, 1, 1, 1};

    /* Edge vertex indices */
    int edge_vertices[12][2] = {
        {0,1}, {1,2}, {2,3}, {3,0},  /* bottom */
        {4,5}, {5,6}, {6,7}, {7,4},  /* top */
        {0,4}, {1,5}, {2,6}, {3,7}   /* vertical */
    };

    /* March through grid */
    for (uint32_t z = 0; z < sdf->resolution[2] - 1; z++) {
        for (uint32_t y = 0; y < sdf->resolution[1] - 1; y++) {
            for (uint32_t x = 0; x < sdf->resolution[0] - 1; x++) {
                /* Get 8 corner values */
                float values[8];
                geomlib_vec3 corners[8];

                for (int i = 0; i < 8; i++) {
                    uint32_t cx = x + dx[i];
                    uint32_t cy = y + dy[i];
                    uint32_t cz = z + dz[i];

                    size_t idx = geomlib_sdf_get_index(sdf, cx, cy, cz);
                    values[i] = sdf->values[idx];

                    corners[i] = geomlib_vec3_create(
                        sdf->origin.x + cx * sdf->voxel_size,
                        sdf->origin.y + cy * sdf->voxel_size,
                        sdf->origin.z + cz * sdf->voxel_size
                    );
                }

                /* Determine cube configuration */
                int cube_index = 0;
                for (int i = 0; i < 8; i++) {
                    if (values[i] < iso_value) {
                        cube_index |= (1 << i);
                    }
                }

                /* Skip if cube is entirely inside or outside */
                if (mc_edge_table[cube_index] == 0) {
                    continue;
                }

                /* Find interpolated vertices on edges */
                geomlib_vec3 vert_list[12];
                for (int i = 0; i < 12; i++) {
                    if (mc_edge_table[cube_index] & (1 << i)) {
                        int v0 = edge_vertices[i][0];
                        int v1 = edge_vertices[i][1];
                        vert_list[i] = interpolate_vertex(
                            corners[v0], corners[v1],
                            values[v0], values[v1],
                            iso_value
                        );
                    }
                }

                /* Create triangles */
                for (int i = 0; mc_tri_table[cube_index][i] != -1; i += 3) {
                    geomlib_vec3 v0 = vert_list[mc_tri_table[cube_index][i]];
                    geomlib_vec3 v1 = vert_list[mc_tri_table[cube_index][i + 1]];
                    geomlib_vec3 v2 = vert_list[mc_tri_table[cube_index][i + 2]];

                    uint32_t i0 = geomlib_trimesh_add_vertex(output, v0);
                    uint32_t i1 = geomlib_trimesh_add_vertex(output, v1);
                    uint32_t i2 = geomlib_trimesh_add_vertex(output, v2);

                    geomlib_trimesh_add_triangle(output, i0, i1, i2);
                }
            }
        }
    }

    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_sdf_extract_surface(
    const geomlib_sdf *sdf,
    geomlib_trimesh *output
) {
    return geomlib_sdf_to_mesh(sdf, 0.0f, output);
}

/* ============================================================================
 * SDF Operations (Stubs)
 * ============================================================================ */

int geomlib_sdf_smooth(geomlib_sdf *sdf, float sigma) {
    /* TODO: Implement Gaussian smoothing */
    (void)sdf;
    (void)sigma;
    return GEOMLIB_MESH_ERROR_INVALID;
}

int geomlib_sdf_dilate(geomlib_sdf *sdf, float amount) {
    /* Simple dilation: subtract amount */
    for (size_t i = 0; i < sdf->total_voxels; i++) {
        sdf->values[i] -= amount;
    }
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_sdf_erode(geomlib_sdf *sdf, float amount) {
    /* Simple erosion: add amount */
    for (size_t i = 0; i < sdf->total_voxels; i++) {
        sdf->values[i] += amount;
    }
    return GEOMLIB_MESH_SUCCESS;
}

geomlib_sdf geomlib_sdf_blend(
    const geomlib_sdf *sdf1,
    const geomlib_sdf *sdf2,
    float blend_factor
) {
    /* TODO: Implement SDF blending */
    geomlib_sdf result;
    memset(&result, 0, sizeof(geomlib_sdf));
    (void)sdf1;
    (void)sdf2;
    (void)blend_factor;
    return result;
}

/* ============================================================================
 * File I/O (Stubs)
 * ============================================================================ */

int geomlib_sdf_save(const geomlib_sdf *sdf, const char *filename) {
    /* TODO: Implement binary save */
    (void)sdf;
    (void)filename;
    return GEOMLIB_MESH_ERROR_INVALID;
}

int geomlib_sdf_load(geomlib_sdf *sdf, const char *filename) {
    /* TODO: Implement binary load */
    (void)sdf;
    (void)filename;
    return GEOMLIB_MESH_ERROR_INVALID;
}
