/*
 * geomlib - Signed Distance Field Implementation
 */

#include "geomlib/geomlib_sdf.h"
#include "geomlib/geomlib_feature.h"
#include "geomlib/geomlib_montecarlo.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

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
    params.sparse = false;
    params.algorithm = GEOMLIB_SDF_ALGORITHM_BRUTE_FORCE;
    return params;
}

/* ============================================================================
 * Sparse brick storage
 * ============================================================================ */

#define GEOMLIB_SDF_BRICK_SHIFT 2u
#define GEOMLIB_SDF_BRICK_DIM (1u << GEOMLIB_SDF_BRICK_SHIFT)
#define GEOMLIB_SDF_BRICK_MASK (GEOMLIB_SDF_BRICK_DIM - 1u)
#define GEOMLIB_SDF_BRICK_VOXELS (GEOMLIB_SDF_BRICK_DIM * GEOMLIB_SDF_BRICK_DIM * \
                                  GEOMLIB_SDF_BRICK_DIM)

struct geomlib_sdf_band {
    uint32_t bricks[3];     /* bricks per axis */
    uint32_t table_size;
    uint32_t brick_count;
    uint32_t brick_capacity;
    uint32_t *offsets;      /* table_size entries: brick slot or UINT32_MAX */
    float *values;          /* brick_count * GEOMLIB_SDF_BRICK_VOXELS */
    uint64_t *valid;        /* brick_count * (GEOMLIB_SDF_BRICK_VOXELS / 64) */
};

static void band_free(struct geomlib_sdf_band *band) {
    if (!band) return;
    free(band->offsets);
    free(band->values);
    free(band->valid);
    free(band);
}

static bool band_alloc(struct geomlib_sdf_band **out, const uint32_t resolution[3]) {
    struct geomlib_sdf_band *band = (struct geomlib_sdf_band*)calloc(1, sizeof(*band));
    if (!band) return false;
    for (int i = 0; i < 3; i++) {
        band->bricks[i] = (resolution[i] + GEOMLIB_SDF_BRICK_MASK) >> GEOMLIB_SDF_BRICK_SHIFT;
        if (!band->bricks[i]) band->bricks[i] = 1;
    }
    uint64_t table = (uint64_t)band->bricks[0] * band->bricks[1] * band->bricks[2];
    if (!table || table > UINT32_MAX) {
        free(band);
        return false;
    }
    band->table_size = (uint32_t)table;
    band->offsets = (uint32_t*)malloc(band->table_size * sizeof(uint32_t));
    if (!band->offsets) {
        free(band);
        return false;
    }
    memset(band->offsets, 0xff, band->table_size * sizeof(uint32_t));
    *out = band;
    return true;
}

static uint32_t band_brick_id(const struct geomlib_sdf_band *band,
                              uint32_t x, uint32_t y, uint32_t z) {
    uint32_t bx = x >> GEOMLIB_SDF_BRICK_SHIFT;
    uint32_t by = y >> GEOMLIB_SDF_BRICK_SHIFT;
    uint32_t bz = z >> GEOMLIB_SDF_BRICK_SHIFT;
    return (bz * band->bricks[1] + by) * band->bricks[0] + bx;
}

static bool band_slot(const struct geomlib_sdf_band *band, uint32_t x, uint32_t y, uint32_t z,
                      uint32_t *slot) {
    uint32_t brick = band_brick_id(band, x, y, z);
    uint32_t offset = band->offsets[brick];
    if (offset == UINT32_MAX) return false;
    uint32_t local = ((z & GEOMLIB_SDF_BRICK_MASK) << (2 * GEOMLIB_SDF_BRICK_SHIFT)) |
                     ((y & GEOMLIB_SDF_BRICK_MASK) << GEOMLIB_SDF_BRICK_SHIFT) |
                     (x & GEOMLIB_SDF_BRICK_MASK);
    if (!(band->valid[(size_t)offset * (GEOMLIB_SDF_BRICK_VOXELS / 64) + (local >> 6)] &
          (UINT64_C(1) << (local & 63)))) return false;
    *slot = offset * GEOMLIB_SDF_BRICK_VOXELS + local;
    return true;
}

static bool band_reserve(struct geomlib_sdf_band *band, uint32_t count) {
    if (count <= band->brick_capacity) return true;
    uint32_t capacity = band->brick_capacity ? band->brick_capacity : 8;
    while (capacity < count) {
        if (capacity > UINT32_MAX / 2) return false;
        capacity *= 2;
    }
#if SIZE_MAX < UINT64_MAX
    if ((uint64_t)capacity * GEOMLIB_SDF_BRICK_VOXELS * sizeof(float) > (uint64_t)SIZE_MAX ||
        (uint64_t)capacity * (GEOMLIB_SDF_BRICK_VOXELS / 64) * sizeof(uint64_t) >
            (uint64_t)SIZE_MAX)
        return false;
#endif
    float *values = (float*)realloc(band->values,
        (size_t)capacity * GEOMLIB_SDF_BRICK_VOXELS * sizeof(float));
    if (!values) return false;
    band->values = values;
    uint64_t *valid = (uint64_t*)realloc(band->valid,
        (size_t)capacity * (GEOMLIB_SDF_BRICK_VOXELS / 64) * sizeof(uint64_t));
    if (!valid) return false;
    band->valid = valid;
    band->brick_capacity = capacity;
    return true;
}

static bool band_write_slot(struct geomlib_sdf_band *band, uint32_t x, uint32_t y, uint32_t z,
                            float value) {
    uint32_t brick = band_brick_id(band, x, y, z);
    if (band->offsets[brick] == UINT32_MAX) {
        if (!band_reserve(band, band->brick_count + 1)) return false;
        uint32_t offset = band->brick_count++;
        band->offsets[brick] = offset;
        memset(band->valid + (size_t)offset * (GEOMLIB_SDF_BRICK_VOXELS / 64), 0,
               (GEOMLIB_SDF_BRICK_VOXELS / 64) * sizeof(uint64_t));
    }
    uint32_t local = ((z & GEOMLIB_SDF_BRICK_MASK) << (2 * GEOMLIB_SDF_BRICK_SHIFT)) |
                     ((y & GEOMLIB_SDF_BRICK_MASK) << GEOMLIB_SDF_BRICK_SHIFT) |
                     (x & GEOMLIB_SDF_BRICK_MASK);
    uint32_t offset = band->offsets[brick];
    band->valid[(size_t)offset * (GEOMLIB_SDF_BRICK_VOXELS / 64) + (local >> 6)] |=
        UINT64_C(1) << (local & 63);
    band->values[(size_t)offset * GEOMLIB_SDF_BRICK_VOXELS + local] = value;
    return true;
}

static bool band_read_slot(const struct geomlib_sdf_band *band, uint32_t x, uint32_t y, uint32_t z,
                           float *value) {
    uint32_t slot;
    if (!band_slot(band, x, y, z, &slot)) return false;
    *value = band->values[slot];
    return true;
}


/* ============================================================================
 * Utilities
 * ============================================================================ */

static bool finite_point(geomlib_vec3 point) {
    return isfinite(point.x) && isfinite(point.y) && isfinite(point.z);
}

static size_t voxel_count(const uint32_t resolution[3]) {
    size_t count = 1;
    for (int i = 0; i < 3; i++) {
        if (!resolution[i] || resolution[i] > INT32_MAX ||
            count > SIZE_MAX / sizeof(float) / resolution[i]) return 0;
        count *= resolution[i];
    }
    return count;
}

static bool valid_grid(const geomlib_sdf *sdf) {
    if (!sdf) return false;
    bool storage = sdf->values != NULL || (sdf->is_sparse && sdf->band != NULL);
    return storage && sdf->total_voxels &&
           sdf->total_voxels == voxel_count(sdf->resolution) &&
           sdf->voxel_size > 0.0f && isfinite(sdf->voxel_size) && finite_point(sdf->origin);
}

static bool sdf_defined(const geomlib_sdf *sdf, uint32_t x, uint32_t y, uint32_t z) {
    if (!sdf || x >= sdf->resolution[0] || y >= sdf->resolution[1] || z >= sdf->resolution[2])
        return false;
    if (sdf->is_sparse) {
        uint32_t slot;
        return band_slot(sdf->band, x, y, z, &slot);
    }
    return sdf->values != NULL;
}

static bool sdf_read(const geomlib_sdf *sdf, uint32_t x, uint32_t y, uint32_t z, float *value) {
    if (!sdf || !value || x >= sdf->resolution[0] || y >= sdf->resolution[1] ||
        z >= sdf->resolution[2]) return false;
    if (sdf->is_sparse) return band_read_slot(sdf->band, x, y, z, value);
    if (!sdf->values) return false;
    *value = sdf->values[geomlib_sdf_get_index(sdf, x, y, z)];
    return true;
}

static bool sdf_write(geomlib_sdf *sdf, uint32_t x, uint32_t y, uint32_t z, float value) {
    if (!sdf || x >= sdf->resolution[0] || y >= sdf->resolution[1] || z >= sdf->resolution[2])
        return false;
    if (sdf->is_sparse) return band_write_slot(sdf->band, x, y, z, value);
    if (!sdf->values) return false;
    sdf->values[geomlib_sdf_get_index(sdf, x, y, z)] = value;
    return true;
}

bool geomlib_sdf_is_defined(const geomlib_sdf *sdf, uint32_t x, uint32_t y, uint32_t z) {
    return sdf_defined(sdf, x, y, z);
}

size_t geomlib_sdf_stored_voxels(const geomlib_sdf *sdf) {
    if (!sdf) return 0;
    if (sdf->is_sparse) return sdf->stored_voxels;
    return sdf->values ? sdf->total_voxels : 0;
}


static bool set_grid_bounds(geomlib_sdf *sdf) {
    const uint32_t *r = sdf->resolution;
    sdf->extent = (geomlib_vec3){r[0] * sdf->voxel_size, r[1] * sdf->voxel_size,
                                 r[2] * sdf->voxel_size};
    sdf->bounds.min = sdf->origin;
    sdf->bounds.max = geomlib_vec3_add(sdf->origin, sdf->extent);
    return finite_point(sdf->extent) && finite_point(sdf->bounds.max);
}

static bool allocate_grid(geomlib_sdf *sdf) {
    sdf->total_voxels = voxel_count(sdf->resolution);
    if (!sdf->total_voxels || !(sdf->voxel_size > 0.0f) || !isfinite(sdf->voxel_size) ||
        !finite_point(sdf->origin) || !set_grid_bounds(sdf)) return false;
    if (sdf->is_sparse) return band_alloc(&sdf->band, sdf->resolution);
    sdf->values = (float*)malloc(sdf->total_voxels * sizeof(float));
    if (sdf->values) sdf->stored_voxels = sdf->total_voxels;
    return sdf->values != NULL;
}

static bool compute_grid_resolution(float voxel_size, geomlib_vec3 extent,
                                    uint32_t resolution[3]) {
    const float values[3] = {extent.x, extent.y, extent.z};
    for (int i = 0; i < 3; i++) {
        double cells = ceil((double)values[i] / voxel_size);
        if (!isfinite(cells) || cells < 0.0 || cells > INT32_MAX) return false;
        resolution[i] = cells < 1.0 ? 1 : (uint32_t)cells;
    }
    return true;
}

float geomlib_compute_optimal_voxel_size(
    const geomlib_trimesh *mesh,
    uint32_t target_resolution
) {
    if (!mesh || !mesh->vertex_count || !mesh->vertices || !target_resolution) return 0.0f;
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
    if (!valid_grid(sdf) || !x || !y || !z) return false;
    double coords[3] = {floor(((double)world_pos.x - sdf->origin.x) / sdf->voxel_size),
                        floor(((double)world_pos.y - sdf->origin.y) / sdf->voxel_size),
                        floor(((double)world_pos.z - sdf->origin.z) / sdf->voxel_size)};
    int32_t *indices[3] = {x, y, z};
    bool inside = true;
    for (int i = 0; i < 3; i++) {
        if (!isfinite(coords[i]) || coords[i] < INT32_MIN || coords[i] > INT32_MAX)
            return false;
    }
    for (int i = 0; i < 3; i++) {
        *indices[i] = (int32_t)coords[i];
        if (coords[i] < 0 || coords[i] >= sdf->resolution[i]) inside = false;
    }
    return inside;
}

/* ============================================================================
 * Distance Computation
 * ============================================================================ */

geomlib_vec3 geomlib_find_closest_point_on_mesh(const geomlib_trimesh *mesh,
                                               const geomlib_bvh *bvh, geomlib_vec3 point,
                                               float *distance) {
    float nearest = FLT_MAX;
    geomlib_vec3 closest = {0};
    if (bvh && geomlib_bvh_closest_point(bvh, point, &closest, &nearest, NULL)) {
        if (distance) *distance = nearest;
        return closest;
    }
    if (mesh && mesh->vertices && mesh->indices && mesh->triangle_count)
        closest = geomlib_closest_point_mesh(point, mesh, &nearest);
    if (distance) *distance = nearest;
    return closest;
}

float geomlib_compute_unsigned_distance(const geomlib_trimesh *mesh, const geomlib_bvh *bvh,
                                        geomlib_vec3 point) {
    float distance;
    geomlib_find_closest_point_on_mesh(mesh, bvh, point, &distance);
    return distance;
}

float geomlib_compute_winding_number(const geomlib_trimesh *mesh, geomlib_vec3 point) {
    if (!mesh || !mesh->vertices || !mesh->indices || !finite_point(point)) return 0.0f;
    double sum = 0.0;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3d directions[3];
        double length[3];
        for (int j = 0; j < 3; j++) {
            uint32_t index = mesh->indices[i * 3 + j];
            if (index >= mesh->vertex_count) return 0.0f;
            geomlib_vec3 p = mesh->vertices[index];
            if (!finite_point(p)) return 0.0f;
            directions[j] = (geomlib_vec3d){(double)p.x - point.x, (double)p.y - point.y,
                                           (double)p.z - point.z};
            length[j] = geomlib_vec3d_length(directions[j]);
        }
        geomlib_vec3d a = directions[0], b = directions[1], c = directions[2];
        double numerator = geomlib_vec3d_dot(a, geomlib_vec3d_cross(b, c));
        if (numerator == 0.0) continue;
        double denominator = length[0] * length[1] * length[2] +
            geomlib_vec3d_dot(a, b) * length[2] + geomlib_vec3d_dot(b, c) * length[0] +
            geomlib_vec3d_dot(c, a) * length[1];
        sum += 2.0 * atan2(numerator, denominator);
    }
    return (float)(sum / (4.0 * M_PI));
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

/* Count crossings by repeatedly tracing the closest hit along a fixed ray. */
static bool point_inside_bvh(const geomlib_bvh *bvh, geomlib_vec3 point) {
    if (!bvh) return false;
    geomlib_vec3 direction = geomlib_vec3_normalize(
        geomlib_vec3_create(0.123456f, 0.789012f, 0.345678f));
    float epsilon = 1e-4f;
    geomlib_vec3 origin = point;
    bool inside = false;
    for (int hits = 0; hits < 256; hits++) {
        geomlib_ray ray = geomlib_ray_create(origin, direction);
        geomlib_ray_hit record;
        if (!geomlib_bvh_trace_ray(bvh, &ray, &record) || !record.hit) break;
        inside = !inside;
        origin = geomlib_vec3_add(record.position, geomlib_vec3_mul(direction, epsilon));
    }
    return inside;
}

bool geomlib_is_point_inside_mesh(
    const geomlib_trimesh *mesh,
    const geomlib_bvh *bvh,
    geomlib_vec3 point
) {
    if (!mesh || !mesh->vertices || !mesh->indices || !finite_point(point)) return false;
    if (bvh) return point_inside_bvh(bvh, point);
    /* Cast ray in arbitrary direction and count intersections */
    geomlib_vec3 ray_dir = geomlib_vec3_create(0.123456f, 0.789012f, 0.345678f);
    ray_dir = geomlib_vec3_normalize(ray_dir);

    bool inside = false;

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 v0 = mesh->vertices[mesh->indices[i * 3 + 0]];
        geomlib_vec3 v1 = mesh->vertices[mesh->indices[i * 3 + 1]];
        geomlib_vec3 v2 = mesh->vertices[mesh->indices[i * 3 + 2]];

        float t;
        if (ray_triangle_intersect(point, ray_dir, v0, v1, v2, &t)) {
            inside = !inside;
        }
    }

    /* Odd number of intersections = inside */
    return inside;
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

/* Region codes stored in the temporary mask. */
enum { REGION_UNDEFINED = 0, REGION_BAND = 1, REGION_OUTSIDE = 2, REGION_INSIDE = 3 };

static void voxel_center(const geomlib_sdf *sdf, uint32_t x, uint32_t y, uint32_t z,
                         geomlib_vec3 *center) {
    float h = sdf->voxel_size;
    center->x = sdf->origin.x + ((float)x + 0.5f) * h;
    center->y = sdf->origin.y + ((float)y + 0.5f) * h;
    center->z = sdf->origin.z + ((float)z + 0.5f) * h;
}

static uint32_t clamp_index(double value, uint32_t end) {
    if (value < 0.0) return 0;
    if (value > (double)end) return end;
    return (uint32_t)value;
}

/* Exact distances (and band membership) for voxel centers near any triangle. */
static void build_band(const geomlib_trimesh *mesh, geomlib_sdf *sdf, float band_radius,
                       uint8_t *mask, float *dist) {
    size_t total = sdf->total_voxels;
    for (size_t i = 0; i < total; i++) {
        dist[i] = FLT_MAX;
        mask[i] = REGION_UNDEFINED;
    }
    float h = sdf->voxel_size;
    const uint32_t rx = sdf->resolution[0] - 1;
    const uint32_t ry = sdf->resolution[1] - 1;
    const uint32_t rz = sdf->resolution[2] - 1;
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        geomlib_vec3 v0 = mesh->vertices[mesh->indices[t * 3 + 0]];
        geomlib_vec3 v1 = mesh->vertices[mesh->indices[t * 3 + 1]];
        geomlib_vec3 v2 = mesh->vertices[mesh->indices[t * 3 + 2]];
        float minx = fminf(v0.x, fminf(v1.x, v2.x)) - band_radius;
        float miny = fminf(v0.y, fminf(v1.y, v2.y)) - band_radius;
        float minz = fminf(v0.z, fminf(v1.z, v2.z)) - band_radius;
        float maxx = fmaxf(v0.x, fmaxf(v1.x, v2.x)) + band_radius;
        float maxy = fmaxf(v0.y, fmaxf(v1.y, v2.y)) + band_radius;
        float maxz = fmaxf(v0.z, fmaxf(v1.z, v2.z)) + band_radius;
        double lo_x = ceil(((double)minx - sdf->origin.x) / h - 0.5);
        double lo_y = ceil(((double)miny - sdf->origin.y) / h - 0.5);
        double lo_z = ceil(((double)minz - sdf->origin.z) / h - 0.5);
        double hi_x = floor(((double)maxx - sdf->origin.x) / h - 0.5);
        double hi_y = floor(((double)maxy - sdf->origin.y) / h - 0.5);
        double hi_z = floor(((double)maxz - sdf->origin.z) / h - 0.5);
        if (hi_x < 0 || hi_y < 0 || hi_z < 0 || lo_x > (double)rx ||
            lo_y > (double)ry || lo_z > (double)rz) continue;
        uint32_t x0 = clamp_index(lo_x, rx), x1 = clamp_index(hi_x, rx);
        uint32_t y0 = clamp_index(lo_y, ry), y1 = clamp_index(hi_y, ry);
        uint32_t z0 = clamp_index(lo_z, rz), z1 = clamp_index(hi_z, rz);
        for (uint32_t z = z0; z <= z1; z++) {
            for (uint32_t y = y0; y <= y1; y++) {
                for (uint32_t x = x0; x <= x1; x++) {
                    geomlib_vec3 center;
                    voxel_center(sdf, x, y, z, &center);
                    geomlib_vec3 closest = geomlib_closest_point_triangle(center, v0, v1, v2);
                    geomlib_vec3 delta = geomlib_vec3_sub(closest, center);
                    float distance = geomlib_vec3_length(delta);
                    if (distance > band_radius) continue;
                    size_t idx = geomlib_sdf_get_index(sdf, x, y, z);
                    if (distance < dist[idx]) {
                        dist[idx] = distance;
                        mask[idx] = REGION_BAND;
                    }
                }
            }
        }
    }
}

/* Resolve the sign of every band sample. Stores signed values into dist. */
static void sign_band(const geomlib_trimesh *mesh, const geomlib_bvh *bvh, geomlib_sdf *sdf,
                      uint8_t *mask, float *dist, bool compute_sign, float limit) {
    for (uint32_t z = 0; z < sdf->resolution[2]; z++) {
        for (uint32_t y = 0; y < sdf->resolution[1]; y++) {
            for (uint32_t x = 0; x < sdf->resolution[0]; x++) {
                size_t idx = geomlib_sdf_get_index(sdf, x, y, z);
                if (mask[idx] != REGION_BAND) continue;
                float value = dist[idx];
                if (compute_sign) {
                    geomlib_vec3 center;
                    voxel_center(sdf, x, y, z, &center);
                    bool inside = bvh ? point_inside_bvh(bvh, center)
                                      : geomlib_is_point_inside_mesh(mesh, NULL, center);
                    if (inside) value = -value;
                }
                if (limit > 0.0f) {
                    if (value > limit) value = limit;
                    if (value < -limit) value = -limit;
                }
                dist[idx] = value;
            }
        }
    }
}

/* Classify non-band voxels by flooding the exterior from the grid boundary. */
static bool classify_regions(const geomlib_sdf *sdf, uint8_t *mask) {
    size_t total = sdf->total_voxels;
    uint32_t *queue = (uint32_t*)malloc(total * sizeof(uint32_t));
    if (!queue) return false;
    size_t head = 0, tail = 0;
    const uint32_t nx = sdf->resolution[0], ny = sdf->resolution[1], nz = sdf->resolution[2];
    for (uint32_t z = 0; z < nz; z++) {
        for (uint32_t y = 0; y < ny; y++) {
            for (uint32_t x = 0; x < nx; x++) {
                if (x != 0 && y != 0 && z != 0 && x != nx - 1 && y != ny - 1 && z != nz - 1)
                    continue;
                size_t idx = geomlib_sdf_get_index(sdf, x, y, z);
                if (mask[idx] == REGION_UNDEFINED) {
                    mask[idx] = REGION_OUTSIDE;
                    queue[tail++] = (uint32_t)idx;
                }
            }
        }
    }
    while (head < tail) {
        uint32_t idx = queue[head++];
        uint32_t x = idx % nx;
        uint32_t y = (idx / nx) % ny;
        uint32_t z = idx / ((size_t)nx * ny);
        const int offsets[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
        for (int i = 0; i < 6; i++) {
            int64_t qx = (int64_t)x + offsets[i][0];
            int64_t qy = (int64_t)y + offsets[i][1];
            int64_t qz = (int64_t)z + offsets[i][2];
            if (qx < 0 || qy < 0 || qz < 0 || qx >= (int64_t)nx || qy >= (int64_t)ny ||
                qz >= (int64_t)nz) continue;
            size_t nidx = geomlib_sdf_get_index(sdf, (uint32_t)qx, (uint32_t)qy, (uint32_t)qz);
            if (mask[nidx] == REGION_UNDEFINED) {
                mask[nidx] = REGION_OUTSIDE;
                queue[tail++] = (uint32_t)nidx;
            }
        }
    }
    free(queue);
    for (size_t i = 0; i < total; i++)
        if (mask[i] == REGION_UNDEFINED) mask[i] = REGION_INSIDE;
    return true;
}

static float update_eikonal(float a, bool has_a, float b, bool has_b, float c, bool has_c,
                            float h) {
    float d[3];
    int n = 0;
    if (has_a) d[n++] = a;
    if (has_b) d[n++] = b;
    if (has_c) d[n++] = c;
    if (!n) return FLT_MAX;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (d[j] < d[i]) { float swap = d[i]; d[i] = d[j]; d[j] = swap; }
    float u = d[0] + h;
    if (n >= 2 && u > d[1]) {
        float delta = 2.0f * h * h - (d[0] - d[1]) * (d[0] - d[1]);
        if (delta < 0.0f) delta = 0.0f;
        u = 0.5f * (d[0] + d[1] + sqrtf(delta));
        if (n >= 3 && u > d[2]) {
            float sum = d[0] + d[1] + d[2];
            float sum2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            float disc = sum * sum - 3.0f * sum2 + 3.0f * h * h;
            if (disc < 0.0f) disc = 0.0f;
            u = (sum + sqrtf(disc)) / 3.0f;
        }
    }
    return u;
}

/* Godunov upwind sweeps on |grad u| = 1 using fixed band seeds. */
static void sweep_unsigned(float *u, const uint8_t *mask, const geomlib_sdf *sdf) {
    const uint32_t nx = sdf->resolution[0], ny = sdf->resolution[1], nz = sdf->resolution[2];
    const float h = sdf->voxel_size;
    for (int round = 0; round < 4; round++) {
        double max_change = 0.0;
        for (int sweep = 0; sweep < 8; sweep++) {
            int sx = (sweep & 1) ? -1 : 1;
            int sy = (sweep & 2) ? -1 : 1;
            int sz = (sweep & 4) ? -1 : 1;
            int x_begin = sx > 0 ? 0 : (int)nx - 1;
            int x_end = sx > 0 ? (int)nx : -1;
            int y_begin = sy > 0 ? 0 : (int)ny - 1;
            int y_end = sy > 0 ? (int)ny : -1;
            int z_begin = sz > 0 ? 0 : (int)nz - 1;
            int z_end = sz > 0 ? (int)nz : -1;
            for (int z = z_begin; z != z_end; z += sz) {
                for (int y = y_begin; y != y_end; y += sy) {
                    for (int x = x_begin; x != x_end; x += sx) {
                        size_t idx = geomlib_sdf_get_index(sdf, (uint32_t)x, (uint32_t)y,
                                                           (uint32_t)z);
                        if (mask[idx] == REGION_BAND) continue;
                        float values[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
                        bool known[3] = {false, false, false};
                        for (int axis = 0; axis < 3; axis++) {
                            for (int side = -1; side <= 1; side += 2) {
                                int q[3] = {x, y, z};
                                q[axis] += side;
                                if (q[axis] < 0 || q[axis] >= (int)sdf->resolution[axis]) continue;
                                size_t nidx = geomlib_sdf_get_index(sdf, (uint32_t)q[0],
                                    (uint32_t)q[1], (uint32_t)q[2]);
                                if (u[nidx] < values[axis]) values[axis] = u[nidx];
                            }
                            known[axis] = values[axis] < FLT_MAX;
                        }
                        float candidate = update_eikonal(values[0], known[0], values[1],
                            known[1], values[2], known[2], h);
                        if (candidate < u[idx]) {
                            double change = (double)u[idx] - candidate;
                            if (change > max_change) max_change = change;
                            u[idx] = candidate;
                        }
                    }
                }
            }
        }
        if (max_change < 1e-4 * (double)h) break;
    }
}

geomlib_sdf geomlib_mesh_to_sdf(
    const geomlib_trimesh *mesh,
    const geomlib_sdf_params *params
) {
    geomlib_sdf sdf;
    memset(&sdf, 0, sizeof(geomlib_sdf));

    if (!mesh || !params || !mesh->vertices || !mesh->indices || !mesh->triangle_count ||
        mesh->triangle_count > UINT32_MAX / 3 || mesh->triangle_count > SIZE_MAX / (3 * sizeof(uint32_t)) ||
        !isfinite(params->voxel_size) || !isfinite(params->padding) || params->padding < 0.0f ||
        (params->narrow_band && (!isfinite(params->band_width) || params->band_width < 0.0f)))
        return sdf;
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        if (!finite_point(mesh->vertices[i])) return sdf;
    }
    for (size_t i = 0; i < mesh->triangle_count * 3; i++) {
        if (mesh->indices[i] >= mesh->vertex_count) return sdf;
    }

    float voxel_size = params->voxel_size;
    if (voxel_size <= 0.0f) voxel_size = geomlib_compute_optimal_voxel_size(mesh, params->resolution);
    if (!(voxel_size > 0.0f) || !isfinite(voxel_size)) return sdf;
    sdf.voxel_size = voxel_size;

    geomlib_aabb bounds = geomlib_trimesh_compute_aabb(mesh);
    float padding = params->padding * voxel_size;
    geomlib_vec3 offset = {padding, padding, padding};
    sdf.origin = geomlib_vec3_sub(bounds.min, offset);
    geomlib_vec3 extent = geomlib_vec3_add(geomlib_vec3_sub(bounds.max, bounds.min),
                                          geomlib_vec3_mul(offset, 2.0f));
    sdf.is_narrow_band = params->narrow_band;
    sdf.narrow_band_width = params->band_width;
    sdf.is_sparse = params->sparse && params->narrow_band;
    if (!compute_grid_resolution(voxel_size, extent, sdf.resolution) || !allocate_grid(&sdf))
        return (geomlib_sdf){0};

    /* Build BVH if requested */
    geomlib_bvh bvh_data;
    geomlib_bvh *bvh = NULL;
    if (params->use_bvh) {
        if (geomlib_bvh_build(&bvh_data, mesh, GEOMLIB_BVH_QUALITY_MEDIUM) == GEOMLIB_MESH_SUCCESS) {
            bvh = &bvh_data;
        }
    }

    bool sweep = params->algorithm == GEOMLIB_SDF_ALGORITHM_FAST_SWEEPING;
    float limit = params->narrow_band ? params->band_width * voxel_size : 0.0f;

    if (!sdf.is_sparse && !sweep) {
        /* Dense brute force. Narrow-band values are clamped. */
        for (uint32_t z = 0; z < sdf.resolution[2]; z++) {
            for (uint32_t y = 0; y < sdf.resolution[1]; y++) {
                for (uint32_t x = 0; x < sdf.resolution[0]; x++) {
                    geomlib_vec3 voxel_center_point;
                    voxel_center(&sdf, x, y, z, &voxel_center_point);
                    float distance;
                    if (params->compute_sign) {
                        distance = geomlib_compute_signed_distance(mesh, bvh, voxel_center_point);
                    } else {
                        distance = geomlib_compute_unsigned_distance(mesh, bvh, voxel_center_point);
                    }
                    if (params->narrow_band) {
                        if (distance > limit) distance = limit;
                        if (distance < -limit) distance = -limit;
                    }
                    sdf.values[geomlib_sdf_get_index(&sdf, x, y, z)] = distance;
                }
            }
        }
        if (bvh) geomlib_bvh_free(bvh);
        return sdf;
    }

    /* Band path: sparse storage and/or eikonal sweeping. */
    float seed_width = params->narrow_band ? params->band_width : 3.0f;
    float band_radius = (seed_width + 1.0f) * voxel_size;
    size_t total = sdf.total_voxels;
    uint8_t *mask = (uint8_t*)malloc(total);
    float *dist = (float*)malloc(total * sizeof(float));
    if (!mask || !dist) {
        free(mask);
        free(dist);
        geomlib_sdf_free(&sdf);
        if (bvh) geomlib_bvh_free(bvh);
        return (geomlib_sdf){0};
    }

    build_band(mesh, &sdf, band_radius, mask, dist);
    sign_band(mesh, bvh, &sdf, mask, dist, params->compute_sign, limit);

    if (sdf.is_sparse) {
        size_t stored = 0;
        for (uint32_t z = 0; z < sdf.resolution[2]; z++) {
            for (uint32_t y = 0; y < sdf.resolution[1]; y++) {
                for (uint32_t x = 0; x < sdf.resolution[0]; x++) {
                    size_t idx = geomlib_sdf_get_index(&sdf, x, y, z);
                    if (mask[idx] != REGION_BAND) continue;
                    if (!sdf_write(&sdf, x, y, z, dist[idx])) {
                        free(mask);
                        free(dist);
                        geomlib_sdf_free(&sdf);
                        if (bvh) geomlib_bvh_free(bvh);
                        return (geomlib_sdf){0};
                    }
                    stored++;
                }
            }
        }
        sdf.stored_voxels = stored;
        free(mask);
        free(dist);
        if (bvh) geomlib_bvh_free(bvh);
        return sdf;
    }

    if (params->compute_sign) {
        if (!classify_regions(&sdf, mask)) {
            free(mask);
            free(dist);
            geomlib_sdf_free(&sdf);
            if (bvh) geomlib_bvh_free(bvh);
            return (geomlib_sdf){0};
        }
    } else {
        for (size_t i = 0; i < total; i++)
            if (mask[i] != REGION_BAND) mask[i] = REGION_OUTSIDE;
    }

    if (!params->narrow_band) {
        float *u = (float*)malloc(total * sizeof(float));
        if (!u) {
            free(mask);
            free(dist);
            geomlib_sdf_free(&sdf);
            if (bvh) geomlib_bvh_free(bvh);
            return (geomlib_sdf){0};
        }
        for (size_t i = 0; i < total; i++)
            u[i] = mask[i] == REGION_BAND ? fabsf(dist[i]) : FLT_MAX;
        sweep_unsigned(u, mask, &sdf);
        for (uint32_t z = 0; z < sdf.resolution[2]; z++) {
            for (uint32_t y = 0; y < sdf.resolution[1]; y++) {
                for (uint32_t x = 0; x < sdf.resolution[0]; x++) {
                    size_t idx = geomlib_sdf_get_index(&sdf, x, y, z);
                    float value;
                    if (mask[idx] == REGION_BAND) {
                        value = dist[idx];
                    } else {
                        value = (mask[idx] == REGION_INSIDE) ? -u[idx] : u[idx];
                    }
                    sdf.values[idx] = value;
                }
            }
        }
        free(u);
    } else {
        /* Clamped narrow band: far samples take the signed band limit. */
        for (uint32_t z = 0; z < sdf.resolution[2]; z++) {
            for (uint32_t y = 0; y < sdf.resolution[1]; y++) {
                for (uint32_t x = 0; x < sdf.resolution[0]; x++) {
                    size_t idx = geomlib_sdf_get_index(&sdf, x, y, z);
                    sdf.values[idx] = mask[idx] == REGION_BAND ? dist[idx]
                        : (mask[idx] == REGION_INSIDE ? -limit : limit);
                }
            }
        }
    }

    free(mask);
    free(dist);
    if (bvh) geomlib_bvh_free(bvh);
    return sdf;
}

void geomlib_sdf_free(geomlib_sdf *sdf) {
    if (!sdf) return;
    free(sdf->values);
    band_free(sdf->band);
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
    float value;
    if (!valid_grid(sdf) || !sdf_read(sdf, x, y, z, &value)) return FLT_MAX;
    return value;
}

static float sample_grid(const geomlib_sdf *sdf, const double coordinate[3], bool clamp) {
    uint32_t lo[3], hi[3];
    double t[3];
    for (int i = 0; i < 3; i++) {
        double end = sdf->resolution[i] - 1;
        double value = coordinate[i];
        if (!isfinite(value)) return FLT_MAX;
        /* Coordinates are computed from float sample positions at the public entry point. */
        double tolerance = 4.0 * FLT_EPSILON * fmax(1.0, end);
        if (!clamp && (value < -tolerance || value > end + tolerance)) return FLT_MAX;
        value = fmax(0.0, fmin(value, end));
        lo[i] = (uint32_t)value;
        hi[i] = lo[i] < end ? lo[i] + 1 : lo[i];
        t[i] = value - lo[i];
    }
    double planes[2];
    for (int z = 0; z < 2; z++) {
        double rows[2];
        for (int y = 0; y < 2; y++) {
            uint32_t j = y ? hi[1] : lo[1], k = z ? hi[2] : lo[2];
            float a, b;
            if (!sdf_read(sdf, lo[0], j, k, &a) || !sdf_read(sdf, hi[0], j, k, &b))
                return FLT_MAX;
            rows[y] = a * (1.0 - t[0]) + b * t[0];
        }
        planes[z] = rows[0] * (1.0 - t[1]) + rows[1] * t[1];
    }
    return (float)(planes[0] * (1.0 - t[2]) + planes[1] * t[2]);
}

float geomlib_sdf_query_distance(const geomlib_sdf *sdf, geomlib_vec3 point) {
    if (!valid_grid(sdf)) return FLT_MAX;
    const double coordinate[3] = {
        ((double)point.x - sdf->origin.x) / sdf->voxel_size - 0.5,
        ((double)point.y - sdf->origin.y) / sdf->voxel_size - 0.5,
        ((double)point.z - sdf->origin.z) / sdf->voxel_size - 0.5
    };
    return sample_grid(sdf, coordinate, false);
}

geomlib_vec3 geomlib_sdf_query_gradient(
    const geomlib_sdf *sdf,
    geomlib_vec3 point
) {
    /* Central differences with one-sided differences at the grid boundary. */
    if (!valid_grid(sdf)) return (geomlib_vec3){0};
    float h = sdf->voxel_size;

    geomlib_vec3 dx = geomlib_vec3_create(h, 0, 0);
    geomlib_vec3 dy = geomlib_vec3_create(0, h, 0);
    geomlib_vec3 dz = geomlib_vec3_create(0, 0, h);

    float p = geomlib_sdf_query_distance(sdf, point);
    if (p == FLT_MAX) return (geomlib_vec3){0};
    const geomlib_vec3 axes[3] = {dx, dy, dz};
    float derivative[3];
    for (int i = 0; i < 3; i++) {
        float low = geomlib_sdf_query_distance(sdf, geomlib_vec3_sub(point, axes[i]));
        float high = geomlib_sdf_query_distance(sdf, geomlib_vec3_add(point, axes[i]));
        if (low == FLT_MAX && high == FLT_MAX) derivative[i] = 0.0f;
        else if (low == FLT_MAX) derivative[i] = (high - p) / h;
        else if (high == FLT_MAX) derivative[i] = (p - low) / h;
        else derivative[i] = (high - low) / (2.0f * h);
    }
    return (geomlib_vec3){derivative[0], derivative[1], derivative[2]};
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

    if (!valid_grid(sdf)) return stats;
    stats.min_distance = FLT_MAX;
    stats.max_distance = -FLT_MAX;
    double sum = 0.0;
    size_t count = 0;

    for (uint32_t z = 0; z < sdf->resolution[2]; z++) {
        for (uint32_t y = 0; y < sdf->resolution[1]; y++) {
            for (uint32_t x = 0; x < sdf->resolution[0]; x++) {
                float d;
                if (!sdf_read(sdf, x, y, z, &d)) continue;
                count++;
                if (d < stats.min_distance) stats.min_distance = d;
                if (d > stats.max_distance) stats.max_distance = d;
                sum += d;
                if (d < 0.0f) {
                    stats.inside_voxels++;
                } else if (d > 0.0f) {
                    stats.outside_voxels++;
                }
                if (fabsf(d) < sdf->voxel_size) stats.surface_voxels++;
            }
        }
    }

    if (!count) {
        stats.min_distance = 0.0f;
        stats.max_distance = 0.0f;
        return stats;
    }
    stats.avg_distance = (float)(sum / (double)count);
    if (sdf->is_sparse && sdf->band) {
        stats.memory_bytes = (size_t)sdf->band->brick_count * GEOMLIB_SDF_BRICK_VOXELS *
                                 sizeof(float) +
                             (size_t)sdf->band->table_size * sizeof(uint32_t);
    } else {
        stats.memory_bytes = sdf->total_voxels * sizeof(float);
    }

    return stats;
}


/* ============================================================================
 * Marching Cubes Implementation
 * ============================================================================ */

/* Each face is wound outward. Opposite faces start at the same transverse
 * corner, so the ambiguous-face rule agrees across neighboring cubes. */
static const uint8_t cube_corners[8][3] = {
    {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
};
static const uint8_t cube_edges[12][3] = {
    {0,1,0}, {1,2,1}, {3,2,0}, {0,3,1}, {4,5,0}, {5,6,1},
    {7,6,0}, {4,7,1}, {0,4,2}, {1,5,2}, {2,6,2}, {3,7,2}
};
static const uint8_t face_corners[6][4] = {
    {0,3,2,1}, {4,5,6,7}, {0,1,5,4}, {3,7,6,2}, {0,4,7,3}, {1,2,6,5}
};
static const uint8_t face_edges[6][4] = {
    {3,2,1,0}, {4,5,6,7}, {0,9,4,8}, {11,6,10,2}, {8,7,11,3}, {1,10,5,9}
};

static void link_face_edges(int8_t next[12], const double values[8], int face, int a, int b) {
    if (values[face_corners[face][a]] < 0.0) {
        int temp = a;
        a = b;
        b = temp;
    }
    next[face_edges[face][a]] = (int8_t)face_edges[face][b];
}

static void build_cube_contours(const double values[8], int8_t next[12]) {
    for (int e = 0; e < 12; e++) next[e] = -1;
    for (int face = 0; face < 6; face++) {
        int crossings[4], count = 0;
        const uint8_t *corners = face_corners[face];
        for (int i = 0; i < 4; i++) {
            if ((values[corners[i]] < 0.0) != (values[corners[(i + 1) % 4]] < 0.0))
                crossings[count++] = i;
        }
        if (count == 2) {
            link_face_edges(next, values, face, crossings[0], crossings[1]);
        } else if (count == 4) {
            /* The sign at the bilinear saddle selects the pairing. */
            double determinant = values[corners[0]] * values[corners[2]] -
                                 values[corners[1]] * values[corners[3]];
            if (determinant >= 0.0) {
                link_face_edges(next, values, face, 0, 1);
                link_face_edges(next, values, face, 2, 3);
            } else {
                link_face_edges(next, values, face, 3, 0);
                link_face_edges(next, values, face, 1, 2);
            }
        }
    }
}

/* Move an extracted vertex onto a reference sharp edge when close enough. */
typedef struct geomlib_sdf_placement {
    const geomlib_trimesh *reference;
    const geomlib_bvh *ref_bvh;
    const geomlib_feature_edge *features;
    size_t feature_count;
    float radius;
    bool snap_to_surface;
} geomlib_sdf_placement;

static geomlib_vec3 closest_point_segment(geomlib_vec3 point, geomlib_vec3 a, geomlib_vec3 b,
                                          float *distance) {
    geomlib_vec3 ab = geomlib_vec3_sub(b, a);
    float length2 = geomlib_vec3_dot(ab, ab);
    float t = length2 > 0.0f ? geomlib_vec3_dot(geomlib_vec3_sub(point, a), ab) / length2 : 0.0f;
    t = fmaxf(0.0f, fminf(1.0f, t));
    geomlib_vec3 closest = geomlib_vec3_add(a, geomlib_vec3_mul(ab, t));
    *distance = geomlib_vec3_length(geomlib_vec3_sub(point, closest));
    return closest;
}

static geomlib_vec3 place_point(const geomlib_sdf_placement *placement, geomlib_vec3 point) {
    if (!placement) return point;
    float best = placement->radius;
    geomlib_vec3 snapped = point;
    bool found = false;
    for (size_t i = 0; i < placement->feature_count; i++) {
        const geomlib_feature_edge *edge = &placement->features[i];
        if (edge->v0 >= placement->reference->vertex_count ||
            edge->v1 >= placement->reference->vertex_count) continue;
        float distance;
        geomlib_vec3 candidate = closest_point_segment(point,
            placement->reference->vertices[edge->v0],
            placement->reference->vertices[edge->v1], &distance);
        if (distance <= best) {
            best = distance;
            snapped = candidate;
            found = true;
        }
    }
    if (found) return snapped;
    if (placement->snap_to_surface && placement->ref_bvh) {
        geomlib_vec3 closest;
        if (geomlib_bvh_closest_point(placement->ref_bvh, point, &closest, NULL, NULL))
            return closest;
    }
    return point;
}

static int contour_vertex(const geomlib_sdf *sdf, geomlib_trimesh *mesh, uint32_t *cache,
                           const uint32_t cell[3], const double values[8], int edge,
                           const geomlib_sdf_placement *placement, uint32_t *vertex) {
    int a = cube_edges[edge][0], b = cube_edges[edge][1];
    int axis = cube_edges[edge][2];
    double t = values[a] / (values[a] - values[b]);
    int corner = t == 1.0 ? b : a;
    size_t sample = geomlib_sdf_get_index(sdf, cell[0] + cube_corners[corner][0],
        cell[1] + cube_corners[corner][1], cell[2] + cube_corners[corner][2]);
    size_t slot = sample * 4 + ((t == 0.0 || t == 1.0) ? 3 : axis);
    if (cache[slot] == UINT32_MAX) {
        double offset[3];
        for (int i = 0; i < 3; i++) offset[i] = cell[i] + cube_corners[a][i] + 0.5;
        offset[axis] += t;
        geomlib_vec3 point = {
            (float)(sdf->origin.x + offset[0] * sdf->voxel_size),
            (float)(sdf->origin.y + offset[1] * sdf->voxel_size),
            (float)(sdf->origin.z + offset[2] * sdf->voxel_size)
        };
        point = place_point(placement, point);
        if (mesh->vertex_count >= UINT32_MAX) return GEOMLIB_MESH_ERROR_MEMORY;
        uint32_t index = (uint32_t)mesh->vertex_count;
        int result = geomlib_trimesh_add_vertex(mesh, point, NULL, NULL);
        if (result != GEOMLIB_MESH_SUCCESS) return result;
        cache[slot] = index;
    }
    *vertex = cache[slot];
    return GEOMLIB_MESH_SUCCESS;
}

static int add_contour_triangle(geomlib_trimesh *mesh, uint32_t a, uint32_t b, uint32_t c) {
    if (a == b || a == c || b == c) return GEOMLIB_MESH_SUCCESS;
    return geomlib_trimesh_add_triangle(mesh, a, b, c);
}

static int emit_cube_contours(const geomlib_sdf *sdf, geomlib_trimesh *mesh, uint32_t *cache,
                               const uint32_t cell[3], const double values[8],
                               const geomlib_sdf_placement *placement) {
    int8_t next[12];
    build_cube_contours(values, next);
    for (int start = 0; start < 12; start++) {
        if (next[start] < 0) continue;
        uint32_t loop[12];
        size_t count = 0;
        int current = start;
        do {
            uint32_t index;
            int result = contour_vertex(sdf, mesh, cache, cell, values, current, placement,
                                        &index);
            if (result != GEOMLIB_MESH_SUCCESS) return result;
            if (!count || loop[count - 1] != index) loop[count++] = index;
            int following = next[current];
            next[current] = -1;
            current = following;
        } while (current != start);
        if (count > 1 && loop[count - 1] == loop[0]) count--;
        if (count < 3) continue;
        if (count == 3) {
            int result = add_contour_triangle(mesh, loop[0], loop[1], loop[2]);
            if (result != GEOMLIB_MESH_SUCCESS) return result;
            continue;
        }

        /* An interior fan introduces no diagonals on shared cube faces. */
        geomlib_vec3 center = {0};
        for (size_t i = 0; i < count; i++) {
            center = geomlib_vec3_add(center, geomlib_vec3_mul(mesh->vertices[loop[i]], 1.0f / count));
        }
        center = place_point(placement, center);
        if (mesh->vertex_count >= UINT32_MAX) return GEOMLIB_MESH_ERROR_MEMORY;
        uint32_t middle = (uint32_t)mesh->vertex_count;
        int result = geomlib_trimesh_add_vertex(mesh, center, NULL, NULL);
        if (result != GEOMLIB_MESH_SUCCESS) return result;
        for (size_t i = 0; i < count; i++) {
            result = add_contour_triangle(mesh, loop[i], loop[(i + 1) % count], middle);
            if (result != GEOMLIB_MESH_SUCCESS) return result;
        }
    }
    return GEOMLIB_MESH_SUCCESS;
}

static int extract_core(const geomlib_sdf *sdf, float iso_value,
                        const geomlib_sdf_placement *placement, geomlib_trimesh *output) {
    if (!output || !valid_grid(sdf) || !isfinite(iso_value)) return GEOMLIB_MESH_ERROR_INVALID;
    geomlib_trimesh mesh = {0};
    int result = geomlib_trimesh_init(&mesh, 0, 0);
    if (result != GEOMLIB_MESH_SUCCESS) return result;
    if (sdf->resolution[0] < 2 || sdf->resolution[1] < 2 || sdf->resolution[2] < 2) {
        *output = mesh;
        return GEOMLIB_MESH_SUCCESS;
    }
    if (sdf->total_voxels > SIZE_MAX / (4 * sizeof(uint32_t))) {
        geomlib_trimesh_free(&mesh);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    size_t bytes = sdf->total_voxels * 4 * sizeof(uint32_t);
    uint32_t *cache = (uint32_t*)malloc(bytes);
    if (!cache) {
        geomlib_trimesh_free(&mesh);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    memset(cache, 0xff, bytes);

    for (uint32_t z = 0; z + 1 < sdf->resolution[2]; z++) {
        for (uint32_t y = 0; y + 1 < sdf->resolution[1]; y++) {
            for (uint32_t x = 0; x + 1 < sdf->resolution[0]; x++) {
                double values[8];
                int inside = 0;
                bool complete = true;
                for (int i = 0; i < 8; i++) {
                    float sample;
                    if (!sdf_read(sdf, x + cube_corners[i][0], y + cube_corners[i][1],
                                  z + cube_corners[i][2], &sample)) {
                        complete = false;
                        break;
                    }
                    values[i] = (double)sample - iso_value;
                    if (!isfinite(values[i])) {
                        result = GEOMLIB_MESH_ERROR_INVALID;
                        goto finish;
                    }
                    inside += values[i] < 0.0;
                }
                if (!complete || !inside || inside == 8) continue;
                const uint32_t cell[3] = {x, y, z};
                result = emit_cube_contours(sdf, &mesh, cache, cell, values, placement);
                if (result != GEOMLIB_MESH_SUCCESS) goto finish;
            }
        }
    }

finish:
    free(cache);
    if (result == GEOMLIB_MESH_SUCCESS) *output = mesh;
    else geomlib_trimesh_free(&mesh);
    return result;
}

int geomlib_sdf_to_mesh(const geomlib_sdf *sdf, float iso_value, geomlib_trimesh *output) {
    return extract_core(sdf, iso_value, NULL, output);
}

int geomlib_sdf_extract_surface_features(
    const geomlib_sdf *sdf,
    const geomlib_trimesh *reference,
    float sharp_angle_deg,
    bool snap_to_surface,
    geomlib_trimesh *output
) {
    if (!output || !valid_grid(sdf) || !isfinite(sharp_angle_deg))
        return GEOMLIB_MESH_ERROR_INVALID;
    if (!reference) return extract_core(sdf, 0.0f, NULL, output);

    geomlib_feature_params defaults = geomlib_feature_params_default();
    defaults.sharp_angle_threshold = sharp_angle_deg;
    defaults.detect_ridges_valleys = false;
    defaults.detect_corners = false;
    defaults.mark_boundaries = false;
    geomlib_features features = geomlib_detect_features(reference, &defaults);

    geomlib_bvh ref_bvh;
    geomlib_bvh *bvh = NULL;
    if (geomlib_bvh_build(&ref_bvh, reference, GEOMLIB_BVH_QUALITY_MEDIUM) ==
        GEOMLIB_MESH_SUCCESS) bvh = &ref_bvh;

    geomlib_sdf_placement placement = {
        reference, bvh, features.edges, features.edge_count, 1.5f * sdf->voxel_size,
        snap_to_surface
    };
    int result = extract_core(sdf, 0.0f, &placement, output);
    if (bvh) geomlib_bvh_free(bvh);
    geomlib_features_free(&features);
    return result;
}

int geomlib_sdf_extract_surface(
    const geomlib_sdf *sdf,
    geomlib_trimesh *output
) {
    return geomlib_sdf_to_mesh(sdf, 0.0f, output);
}

/* ============================================================================
 * SDF Operations
 * ============================================================================ */

/* Defined-neighbour Gaussian blur. Undefined samples are ignored and left unset. */
int geomlib_sdf_smooth(geomlib_sdf *sdf, float sigma) {
    if (!valid_grid(sdf) || !(sigma > 0.0f) || !isfinite(sigma))
        return GEOMLIB_MESH_ERROR_INVALID;

    uint32_t longest = sdf->resolution[0];
    for (int i = 1; i < 3; i++) {
        if (sdf->resolution[i] > longest) longest = sdf->resolution[i];
    }
    int64_t radius = (int64_t)fmin(ceil(3.0 * sigma), longest - 1);
    size_t kernel_size = (size_t)(2 * radius + 1);
    if (kernel_size > SIZE_MAX / sizeof(float)) return GEOMLIB_MESH_ERROR_MEMORY;
    float *kernel = (float*)malloc(kernel_size * sizeof(float));
    float *smoothed = (float*)malloc(sdf->total_voxels * sizeof(float));
    float *temp = (float*)malloc(sdf->total_voxels * sizeof(float));
    if (!kernel || !smoothed || !temp) {
        free(kernel);
        free(smoothed);
        free(temp);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    memset(temp, 0, sdf->total_voxels * sizeof(float));
    memset(smoothed, 0, sdf->total_voxels * sizeof(float));
    for (int64_t i = -radius; i <= radius; i++) {
        double offset = (double)i / sigma;
        kernel[i + radius] = (float)exp(-0.5 * offset * offset);
    }

    const uint32_t r[3] = {sdf->resolution[0], sdf->resolution[1], sdf->resolution[2]};
    const float *sources[3] = {NULL, temp, smoothed};
    float *targets[3] = {temp, smoothed, temp};
    for (int axis = 0; axis < 3; axis++) {
        for (uint32_t z = 0; z < r[2]; z++) {
            for (uint32_t y = 0; y < r[1]; y++) {
                for (uint32_t x = 0; x < r[0]; x++) {
                    if (!sdf_defined(sdf, x, y, z)) continue;
                    float sum = 0.0f;
                    float weight_sum = 0.0f;
                    for (int64_t k = -radius; k <= radius; k++) {
                        uint32_t q[3] = {x, y, z};
                        int64_t coordinate = (int64_t)q[axis] + k;
                        if (coordinate < 0 || coordinate >= (int64_t)r[axis]) continue;
                        q[axis] = (uint32_t)coordinate;
                        if (!sdf_defined(sdf, q[0], q[1], q[2])) continue;
                        float value;
                        if (axis == 0) {
                            if (!sdf_read(sdf, q[0], q[1], q[2], &value)) continue;
                        } else {
                            value = sources[axis][geomlib_sdf_get_index(sdf, q[0], q[1], q[2])];
                        }
                        float w = kernel[k + radius];
                        sum += value * w;
                        weight_sum += w;
                    }
                    if (weight_sum > 0.0f)
                        targets[axis][geomlib_sdf_get_index(sdf, x, y, z)] = sum / weight_sum;
                }
            }
        }
    }

    for (uint32_t z = 0; z < r[2]; z++) {
        for (uint32_t y = 0; y < r[1]; y++) {
            for (uint32_t x = 0; x < r[0]; x++) {
                if (!sdf_defined(sdf, x, y, z)) continue;
                sdf_write(sdf, x, y, z, temp[geomlib_sdf_get_index(sdf, x, y, z)]);
            }
        }
    }

    free(kernel);
    free(temp);
    free(smoothed);

    return GEOMLIB_MESH_SUCCESS;
}

static int shift_defined(geomlib_sdf *sdf, float amount) {
    if (!valid_grid(sdf) || !isfinite(amount)) return GEOMLIB_MESH_ERROR_INVALID;
    for (uint32_t z = 0; z < sdf->resolution[2]; z++) {
        for (uint32_t y = 0; y < sdf->resolution[1]; y++) {
            for (uint32_t x = 0; x < sdf->resolution[0]; x++) {
                float value;
                if (!sdf_read(sdf, x, y, z, &value)) continue;
                sdf_write(sdf, x, y, z, value + amount);
            }
        }
    }
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_sdf_dilate(geomlib_sdf *sdf, float amount) {
    return shift_defined(sdf, -amount);
}

int geomlib_sdf_erode(geomlib_sdf *sdf, float amount) {
    return shift_defined(sdf, amount);
}

geomlib_sdf geomlib_sdf_blend(
    const geomlib_sdf *sdf1,
    const geomlib_sdf *sdf2,
    float blend_factor
) {
    geomlib_sdf result;
    memset(&result, 0, sizeof(geomlib_sdf));

    /* Validate inputs */
    /* Sparse grids are not blendable; materialize first with geomlib_sdf_make_dense. */
    if (!valid_grid(sdf1) || !valid_grid(sdf2) || !isfinite(blend_factor) ||
        sdf1->is_sparse || sdf2->is_sparse) {
        return result;
    }

    /* Check that SDFs have compatible dimensions */
    if (sdf1->resolution[0] != sdf2->resolution[0] ||
        sdf1->resolution[1] != sdf2->resolution[1] ||
        sdf1->resolution[2] != sdf2->resolution[2] ||
        fabsf(sdf1->voxel_size - sdf2->voxel_size) > 1e-6f) {
        return result;  /* Incompatible SDFs */
    }

    /* Clamp blend factor to [0, 1] */
    blend_factor = fmaxf(0.0f, fminf(1.0f, blend_factor));

    /* Allocate result SDF */
    result.resolution[0] = sdf1->resolution[0];
    result.resolution[1] = sdf1->resolution[1];
    result.resolution[2] = sdf1->resolution[2];
    result.total_voxels = sdf1->total_voxels;
    result.voxel_size = sdf1->voxel_size;
    result.origin = sdf1->origin;
    result.extent = sdf1->extent;
    result.bounds = sdf1->bounds;
    result.is_narrow_band = sdf1->is_narrow_band || sdf2->is_narrow_band;
    result.narrow_band_width = fmaxf(sdf1->narrow_band_width, sdf2->narrow_band_width);
    result.stored_voxels = result.total_voxels;

    result.values = (float*)malloc(result.total_voxels * sizeof(float));
    if (!result.values) {
        memset(&result, 0, sizeof(geomlib_sdf));
        return result;
    }

    /* Blend distance values: linear interpolation */
    for (size_t i = 0; i < result.total_voxels; i++) {
        result.values[i] = sdf1->values[i] * (1.0f - blend_factor) +
                          sdf2->values[i] * blend_factor;
    }

    return result;
}

int geomlib_sdf_make_dense(geomlib_sdf *sdf) {
    if (!valid_grid(sdf)) return GEOMLIB_MESH_ERROR_INVALID;
    if (!sdf->is_sparse) {
        sdf->stored_voxels = sdf->total_voxels;
        return GEOMLIB_MESH_SUCCESS;
    }
    size_t total = sdf->total_voxels;
    float *dense = (float*)malloc(total * sizeof(float));
    uint8_t *known = (uint8_t*)calloc(total, 1);
    uint32_t *queue = (uint32_t*)malloc(total * sizeof(uint32_t));
    if (!dense || !known || !queue) {
        free(dense);
        free(known);
        free(queue);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    float limit = sdf->narrow_band_width > 0.0f ? sdf->narrow_band_width * sdf->voxel_size : 0.0f;
    size_t head = 0, tail = 0;
    for (uint32_t z = 0; z < sdf->resolution[2]; z++) {
        for (uint32_t y = 0; y < sdf->resolution[1]; y++) {
            for (uint32_t x = 0; x < sdf->resolution[0]; x++) {
                float value;
                size_t idx = geomlib_sdf_get_index(sdf, x, y, z);
                if (!sdf_read(sdf, x, y, z, &value)) {
                    dense[idx] = 0.0f;
                    continue;
                }
                dense[idx] = value;
                known[idx] = 1;
                if (fabsf(value) > limit) limit = fabsf(value);
                queue[tail++] = (uint32_t)idx;
            }
        }
    }
    const uint32_t nx = sdf->resolution[0], ny = sdf->resolution[1], nz = sdf->resolution[2];
    const int offsets[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
    while (head < tail) {
        uint32_t idx = queue[head++];
        uint32_t x = idx % nx;
        uint32_t y = (idx / nx) % ny;
        uint32_t z = idx / ((size_t)nx * ny);
        float sign = dense[idx] < 0.0f ? -1.0f : 1.0f;
        for (int i = 0; i < 6; i++) {
            int64_t qx = (int64_t)x + offsets[i][0];
            int64_t qy = (int64_t)y + offsets[i][1];
            int64_t qz = (int64_t)z + offsets[i][2];
            if (qx < 0 || qy < 0 || qz < 0 || qx >= (int64_t)nx || qy >= (int64_t)ny ||
                qz >= (int64_t)nz) continue;
            size_t nidx = geomlib_sdf_get_index(sdf, (uint32_t)qx, (uint32_t)qy, (uint32_t)qz);
            if (known[nidx]) continue;
            known[nidx] = 1;
            dense[nidx] = sign * limit;
            queue[tail++] = (uint32_t)nidx;
        }
    }
    band_free(sdf->band);
    sdf->band = NULL;
    sdf->is_sparse = false;
    sdf->values = dense;
    sdf->stored_voxels = total;
    free(known);
    free(queue);
    return GEOMLIB_MESH_SUCCESS;
}

geomlib_sdf geomlib_sdf_resize(const geomlib_sdf *sdf, uint32_t new_resolution) {
    geomlib_sdf result = {0};
    if (!valid_grid(sdf) || sdf->is_sparse || !new_resolution || new_resolution > INT32_MAX)
        return result;
    uint32_t largest = sdf->resolution[0];
    for (int i = 1; i < 3; i++) {
        if (sdf->resolution[i] > largest) largest = sdf->resolution[i];
    }
    result.origin = sdf->origin;
    result.voxel_size = (float)((double)sdf->voxel_size * largest / new_resolution);
    for (int i = 0; i < 3; i++) {
        uint64_t numerator = (uint64_t)sdf->resolution[i] * new_resolution;
        result.resolution[i] = (uint32_t)((numerator + largest - 1) / largest);
    }
    result.is_narrow_band = sdf->is_narrow_band;
    result.narrow_band_width = (float)((double)sdf->narrow_band_width * new_resolution / largest);
    if (!allocate_grid(&result)) return (geomlib_sdf){0};
    result.stored_voxels = result.total_voxels;

    double scale = (double)largest / new_resolution;
    for (uint32_t z = 0; z < result.resolution[2]; z++) {
        for (uint32_t y = 0; y < result.resolution[1]; y++) {
            for (uint32_t x = 0; x < result.resolution[0]; x++) {
                const double coordinate[3] = {(x + 0.5) * scale - 0.5, (y + 0.5) * scale - 0.5,
                                              (z + 0.5) * scale - 0.5};
                result.values[geomlib_sdf_get_index(&result, x, y, z)] = sample_grid(sdf, coordinate, true);
            }
        }
    }
    return result;
}

geomlib_sdf geomlib_sdf_extract_region(const geomlib_sdf *sdf, geomlib_aabb region) {
    geomlib_sdf result = {0};
    if (!valid_grid(sdf) || !finite_point(region.min) || !finite_point(region.max)) return result;
    const float min[3] = {region.min.x, region.min.y, region.min.z};
    const float max[3] = {region.max.x, region.max.y, region.max.z};
    const float origin[3] = {sdf->origin.x, sdf->origin.y, sdf->origin.z};
    uint32_t start[3];
    for (int i = 0; i < 3; i++) {
        double lo = ceil(((double)min[i] - origin[i]) / sdf->voxel_size - 0.5);
        double hi = floor(((double)max[i] - origin[i]) / sdf->voxel_size - 0.5);
        lo = fmax(0, lo);
        hi = fmin(sdf->resolution[i] - 1, hi);
        if (lo > hi) return (geomlib_sdf){0};
        start[i] = (uint32_t)lo;
        result.resolution[i] = (uint32_t)hi - start[i] + 1;
    }
    result.origin = (geomlib_vec3){origin[0] + start[0] * sdf->voxel_size,
                                  origin[1] + start[1] * sdf->voxel_size,
                                  origin[2] + start[2] * sdf->voxel_size};
    result.voxel_size = sdf->voxel_size;
    result.is_narrow_band = sdf->is_narrow_band;
    result.narrow_band_width = sdf->narrow_band_width;
    result.is_sparse = sdf->is_sparse;
    if (!allocate_grid(&result)) return (geomlib_sdf){0};
    size_t stored = 0;
    for (uint32_t z = 0; z < result.resolution[2]; z++) {
        for (uint32_t y = 0; y < result.resolution[1]; y++) {
            for (uint32_t x = 0; x < result.resolution[0]; x++) {
                float value;
                if (!sdf_read(sdf, start[0] + x, start[1] + y, start[2] + z, &value)) continue;
                if (!sdf_write(&result, x, y, z, value)) {
                    geomlib_sdf_free(&result);
                    return (geomlib_sdf){0};
                }
                stored++;
            }
        }
    }
    result.stored_voxels = result.is_sparse ? stored : result.total_voxels;
    if (result.is_sparse && !stored) {
        geomlib_sdf_free(&result);
        return (geomlib_sdf){0};
    }
    return result;
}

/* ============================================================================
 * File I/O
 * ============================================================================ */

int geomlib_sdf_save(const geomlib_sdf *sdf, const char *filename) {
    if (!valid_grid(sdf) || !filename) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Dense grids keep the SDF1 layout; sparse bands use SDF2. */
    const char magic[4] = {'S', 'D', 'F', sdf->is_sparse ? '2' : '1'};
    if (fwrite(magic, 1, 4, fp) != 4 ||
        fwrite(&sdf->resolution[0], sizeof(uint32_t), 3, fp) != 3 ||
        fwrite(&sdf->voxel_size, sizeof(float), 1, fp) != 1 ||
        fwrite(&sdf->origin, sizeof(geomlib_vec3), 1, fp) != 1 ||
        fwrite(&sdf->extent, sizeof(geomlib_vec3), 1, fp) != 1) {
        fclose(fp);
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    uint32_t is_narrow_band = sdf->is_narrow_band ? 1 : 0;
    if (fwrite(&is_narrow_band, sizeof(uint32_t), 1, fp) != 1 ||
        fwrite(&sdf->narrow_band_width, sizeof(float), 1, fp) != 1) {
        fclose(fp);
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    if (sdf->is_sparse) {
        const struct geomlib_sdf_band *band = sdf->band;
        uint32_t brick_shift = GEOMLIB_SDF_BRICK_SHIFT;
        if (fwrite(&brick_shift, sizeof(uint32_t), 1, fp) != 1 ||
            fwrite(&band->brick_count, sizeof(uint32_t), 1, fp) != 1 ||
            fwrite(&band->table_size, sizeof(uint32_t), 1, fp) != 1 ||
            (band->table_size &&
             fwrite(band->offsets, sizeof(uint32_t), band->table_size, fp) != band->table_size) ||
            (band->brick_count &&
             fwrite(band->valid, sizeof(uint64_t),
                    (size_t)band->brick_count * (GEOMLIB_SDF_BRICK_VOXELS / 64), fp) !=
                 (size_t)band->brick_count * (GEOMLIB_SDF_BRICK_VOXELS / 64)) ||
            (band->brick_count &&
             fwrite(band->values, sizeof(float),
                    (size_t)band->brick_count * GEOMLIB_SDF_BRICK_VOXELS, fp) !=
                 (size_t)band->brick_count * GEOMLIB_SDF_BRICK_VOXELS)) {
            fclose(fp);
            return GEOMLIB_MESH_ERROR_INVALID;
        }
    } else if (fwrite(sdf->values, sizeof(float), sdf->total_voxels, fp) != sdf->total_voxels) {
        fclose(fp);
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    fclose(fp);
    return GEOMLIB_MESH_SUCCESS;
}

static bool load_sdf1(FILE *fp, geomlib_sdf *loaded) {
    uint32_t is_narrow_band;
    if (fread(loaded->resolution, sizeof(uint32_t), 3, fp) != 3 ||
        fread(&loaded->voxel_size, sizeof(float), 1, fp) != 1 ||
        fread(&loaded->origin, sizeof(geomlib_vec3), 1, fp) != 1 ||
        fread(&loaded->extent, sizeof(geomlib_vec3), 1, fp) != 1 ||
        fread(&is_narrow_band, sizeof(uint32_t), 1, fp) != 1 ||
        fread(&loaded->narrow_band_width, sizeof(float), 1, fp) != 1 ||
        !finite_point(loaded->extent) || !isfinite(loaded->narrow_band_width)) {
        return false;
    }
    loaded->is_narrow_band = is_narrow_band != 0;
    geomlib_vec3 saved_extent = loaded->extent;
    if (!allocate_grid(loaded)) return false;
    loaded->extent = saved_extent;
    loaded->stored_voxels = loaded->total_voxels;
    if (fread(loaded->values, sizeof(float), loaded->total_voxels, fp) != loaded->total_voxels)
        return false;
    return true;
}

static bool load_sdf2(FILE *fp, geomlib_sdf *loaded) {
    uint32_t is_narrow_band, brick_shift, brick_count, table_size;
    if (fread(loaded->resolution, sizeof(uint32_t), 3, fp) != 3 ||
        fread(&loaded->voxel_size, sizeof(float), 1, fp) != 1 ||
        fread(&loaded->origin, sizeof(geomlib_vec3), 1, fp) != 1 ||
        fread(&loaded->extent, sizeof(geomlib_vec3), 1, fp) != 1 ||
        fread(&is_narrow_band, sizeof(uint32_t), 1, fp) != 1 ||
        fread(&loaded->narrow_band_width, sizeof(float), 1, fp) != 1 ||
        fread(&brick_shift, sizeof(uint32_t), 1, fp) != 1 ||
        fread(&brick_count, sizeof(uint32_t), 1, fp) != 1 ||
        fread(&table_size, sizeof(uint32_t), 1, fp) != 1)
        return false;
    if (brick_shift != GEOMLIB_SDF_BRICK_SHIFT || !finite_point(loaded->extent) ||
        !isfinite(loaded->narrow_band_width)) return false;
    uint32_t expected_bricks[3];
    uint64_t expected_table = 1;
    for (int i = 0; i < 3; i++) {
        expected_bricks[i] = (loaded->resolution[i] + GEOMLIB_SDF_BRICK_MASK) >>
                             GEOMLIB_SDF_BRICK_SHIFT;
        if (!expected_bricks[i]) expected_bricks[i] = 1;
        expected_table *= expected_bricks[i];
    }
    if (!expected_table || expected_table > UINT32_MAX || table_size != (uint32_t)expected_table ||
        !brick_count || brick_count > table_size)
        return false;
    loaded->total_voxels = voxel_count(loaded->resolution);
    if (!loaded->total_voxels) return false;
    geomlib_vec3 saved_extent = loaded->extent;
    if (!set_grid_bounds(loaded)) return false;
    loaded->extent = saved_extent;

    struct geomlib_sdf_band *band = (struct geomlib_sdf_band*)calloc(1, sizeof(*band));
    if (!band) return false;
    for (int i = 0; i < 3; i++) band->bricks[i] = expected_bricks[i];
    band->table_size = table_size;
    band->brick_count = brick_count;
    band->brick_capacity = brick_count;
    band->offsets = (uint32_t*)malloc((size_t)table_size * sizeof(uint32_t));
    band->values = (float*)malloc((size_t)brick_count * GEOMLIB_SDF_BRICK_VOXELS * sizeof(float));
    band->valid = (uint64_t*)malloc((size_t)brick_count *
                                    (GEOMLIB_SDF_BRICK_VOXELS / 64) * sizeof(uint64_t));
    bool ok = band->offsets && band->values && band->valid;
    if (ok) {
        ok = fread(band->offsets, sizeof(uint32_t), table_size, fp) == table_size &&
             fread(band->valid, sizeof(uint64_t),
                   (size_t)brick_count * (GEOMLIB_SDF_BRICK_VOXELS / 64), fp) ==
                 (size_t)brick_count * (GEOMLIB_SDF_BRICK_VOXELS / 64) &&
             fread(band->values, sizeof(float),
                   (size_t)brick_count * GEOMLIB_SDF_BRICK_VOXELS, fp) ==
                 (size_t)brick_count * GEOMLIB_SDF_BRICK_VOXELS;
    }
    if (!ok) {
        band_free(band);
        return false;
    }
    for (uint32_t i = 0; i < table_size; i++) {
        if (band->offsets[i] != UINT32_MAX && band->offsets[i] >= brick_count) {
            band_free(band);
            return false;
        }
    }
    size_t stored = 0;
    for (size_t i = 0; i < (size_t)brick_count * (GEOMLIB_SDF_BRICK_VOXELS / 64); i++) {
        uint64_t bits = band->valid[i];
        while (bits) { stored += (size_t)(bits & 1u); bits >>= 1; }
    }
    loaded->is_sparse = true;
    loaded->is_narrow_band = is_narrow_band != 0;
    loaded->band = band;
    loaded->stored_voxels = stored;
    return true;
}

int geomlib_sdf_load(geomlib_sdf *sdf, const char *filename) {
    if (!sdf || !filename) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    char magic[4];
    if (fread(magic, 1, 4, fp) != 4 ||
        magic[0] != 'S' || magic[1] != 'D' || magic[2] != 'F' ||
        (magic[3] != '1' && magic[3] != '2')) {
        fclose(fp);
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    geomlib_sdf loaded = {0};
    bool ok = magic[3] == '1' ? load_sdf1(fp, &loaded) : load_sdf2(fp, &loaded);
    fclose(fp);
    if (!ok) {
        geomlib_sdf_free(&loaded);
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    *sdf = loaded;
    return GEOMLIB_MESH_SUCCESS;
}
