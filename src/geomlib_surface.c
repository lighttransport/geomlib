/*
 * geomlib - Surface Reconstruction Implementation
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 */

#include "geomlib/geomlib_surface.h"
#include "geomlib/geomlib_random.h"
#include "geomlib/geomlib_sdf.h"
#include "geomlib/geomlib_tetramesh.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <stdint.h>

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

    if (!initial_capacity) initial_capacity = 16;
    if (initial_capacity > SIZE_MAX / sizeof(geomlib_vec3)) return GEOMLIB_MESH_ERROR_INVALID;
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
        if (pc->point_capacity > SIZE_MAX / (2 * sizeof(geomlib_vec3))) return (size_t)-1;
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

    if (pc->normals) {
        pc->normals[idx] = normal ? *normal : (geomlib_vec3){0};
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

static geomlib_vec3d point_as_double(geomlib_vec3 point) {
    return (geomlib_vec3d){point.x, point.y, point.z};
}

static geomlib_vec3 point_as_float(geomlib_vec3d point) {
    return (geomlib_vec3){(float)point.x, (float)point.y, (float)point.z};
}

static bool finite_point(geomlib_vec3 p) {
    return isfinite(p.x) && isfinite(p.y) && isfinite(p.z);
}

typedef struct voxel_sample {
    int64_t cell[3];
    size_t index;
} voxel_sample;

static int compare_voxels(const void *left, const void *right) {
    const voxel_sample *a = left, *b = right;
    for (int i = 0; i < 3; i++) {
        if (a->cell[i] != b->cell[i]) return a->cell[i] < b->cell[i] ? -1 : 1;
    }
    return 0;
}

static int compare_voxel_samples(const void *left, const void *right) {
    int order = compare_voxels(left, right);
    if (order) return order;
    const voxel_sample *a = left, *b = right;
    return (a->index > b->index) - (a->index < b->index);
}

int geomlib_pointcloud_downsample(const geomlib_pointcloud *pc, float voxel_size,
                                  geomlib_pointcloud *output) {
    if (!pc || !output || pc == output || !(voxel_size > 0) || !isfinite(voxel_size) ||
        (pc->point_count && !pc->points) || pc->point_count > SIZE_MAX / sizeof(voxel_sample))
        return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_pointcloud result = {0};
    if (!pc->point_count) {
        int status = geomlib_pointcloud_init(&result, 0);
        if (status == GEOMLIB_MESH_SUCCESS) *output = result;
        return status;
    }
    for (size_t i = 0; i < pc->point_count; i++) {
        if (!finite_point(pc->points[i]) || (pc->normals && !finite_point(pc->normals[i])) ||
            (pc->confidences && !isfinite(pc->confidences[i]))) return GEOMLIB_MESH_ERROR_INVALID;
    }

    voxel_sample *samples = malloc(pc->point_count * sizeof(*samples));
    if (!samples) return GEOMLIB_MESH_ERROR_MEMORY;
    geomlib_vec3 origin = geomlib_pointcloud_bounds(pc).min;
    int status = GEOMLIB_MESH_ERROR_INVALID;
    for (size_t i = 0; i < pc->point_count; i++) {
        geomlib_vec3 p = pc->points[i];
        double cell[3] = {
            floor(((double)p.x - origin.x) / voxel_size),
            floor(((double)p.y - origin.y) / voxel_size),
            floor(((double)p.z - origin.z) / voxel_size)
        };
        for (int axis = 0; axis < 3; axis++) {
            if (cell[axis] >= 0x1p63) goto finish;
            samples[i].cell[axis] = (int64_t)cell[axis];
        }
        samples[i].index = i;
    }
    qsort(samples, pc->point_count, sizeof(*samples), compare_voxel_samples);
    size_t voxels = 1;
    for (size_t i = 1; i < pc->point_count; i++)
        voxels += compare_voxels(&samples[i - 1], &samples[i]) != 0;

    status = geomlib_pointcloud_init(&result, voxels);
    if (status != GEOMLIB_MESH_SUCCESS) goto finish;
    status = GEOMLIB_MESH_ERROR_MEMORY;
    if (pc->normals) {
        result.normals = malloc(voxels * sizeof(geomlib_vec3));
        if (!result.normals) goto finish;
    }
    if (pc->confidences) {
        result.confidences = malloc(voxels * sizeof(float));
        if (!result.confidences) goto finish;
    }
    for (size_t start = 0; start < pc->point_count;) {
        geomlib_vec3d sum = {0}, normals = {0};
        double confidence = 0;
        size_t end = start;
        do {
            size_t index = samples[end++].index;
            sum = geomlib_vec3d_add(sum, point_as_double(pc->points[index]));
            if (pc->normals) normals = geomlib_vec3d_add(normals, point_as_double(pc->normals[index]));
            if (pc->confidences) confidence += pc->confidences[index];
        } while (end < pc->point_count && compare_voxels(&samples[start], &samples[end]) == 0);

        size_t index = result.point_count++;
        double count = (double)(end - start);
        result.points[index] = point_as_float(geomlib_vec3d_mul(sum, 1 / count));
        double length = geomlib_vec3d_length(normals);
        if (result.normals) result.normals[index] = length > 0
            ? point_as_float(geomlib_vec3d_mul(normals, 1 / length)) : (geomlib_vec3){0};
        if (result.confidences) result.confidences[index] = (float)(confidence / count);
        start = end;
    }
    status = GEOMLIB_MESH_SUCCESS;

finish:
    free(samples);
    if (status == GEOMLIB_MESH_SUCCESS) *output = result;
    else geomlib_pointcloud_free(&result);
    return status;
}

typedef struct normal_neighbor {
    uint32_t index;
    double distance;
} normal_neighbor;

static geomlib_vec3 estimate_pca_normal(const geomlib_pointcloud *pc, size_t point,
                                        const normal_neighbor *neighbors, size_t count) {
    geomlib_vec3d origin = point_as_double(pc->points[point]);
    geomlib_vec3d centroid = origin;
    for (size_t i = 0; i < count; i++)
        centroid = geomlib_vec3d_add(centroid, point_as_double(pc->points[neighbors[i].index]));
    centroid = geomlib_vec3d_mul(centroid, 1.0 / (count + 1));

    double covariance[3][3] = {{0}};
    for (size_t i = 0; i <= count; i++) {
        geomlib_vec3d p = i == count ? origin : point_as_double(pc->points[neighbors[i].index]);
        geomlib_vec3d offset = geomlib_vec3d_sub(p, centroid);
        const double v[] = {offset.x, offset.y, offset.z};
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) covariance[j][k] += v[j] * v[k];
        }
    }

    double basis[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    /* Diagonalize the symmetric covariance matrix with Jacobi rotations. */
    for (int sweep = 0; sweep < 32; sweep++) {
        int p = 0, q = 1;
        double scale = 0;
        for (int i = 0; i < 3; i++) {
            scale = fmax(scale, fabs(covariance[i][i]));
            for (int j = i + 1; j < 3; j++) {
                if (fabs(covariance[i][j]) > fabs(covariance[p][q])) {
                    p = i;
                    q = j;
                }
            }
        }
        double off_diagonal = covariance[p][q];
        if (fabs(off_diagonal) <= DBL_EPSILON * scale) break;
        double tau = (covariance[q][q] - covariance[p][p]) / (2 * off_diagonal);
        double tangent = copysign(1.0, tau) / (fabs(tau) + hypot(1.0, tau));
        double c = 1.0 / hypot(1.0, tangent), s = tangent * c;
        covariance[p][p] -= tangent * off_diagonal;
        covariance[q][q] += tangent * off_diagonal;
        covariance[p][q] = covariance[q][p] = 0;
        for (int k = 0; k < 3; k++) {
            if (k != p && k != q) {
                double x = covariance[k][p], y = covariance[k][q];
                covariance[k][p] = covariance[p][k] = c * x - s * y;
                covariance[k][q] = covariance[q][k] = s * x + c * y;
            }
            double x = basis[k][p], y = basis[k][q];
            basis[k][p] = c * x - s * y;
            basis[k][q] = s * x + c * y;
        }
    }

    int smallest = 0, largest = 0;
    for (int i = 1; i < 3; i++) {
        if (covariance[i][i] < covariance[smallest][smallest]) smallest = i;
        if (covariance[i][i] >= covariance[largest][largest]) largest = i;
    }
    int middle = 3 - smallest - largest;
    if (covariance[middle][middle] <= 64 * DBL_EPSILON * covariance[largest][largest])
        return (geomlib_vec3){0};

    double sign = 1, magnitude = 0;
    for (int i = 0; i < 3; i++) {
        if (fabs(basis[i][smallest]) <= magnitude) continue;
        magnitude = fabs(basis[i][smallest]);
        sign = copysign(1, basis[i][smallest]);
    }
    return (geomlib_vec3){(float)(sign * basis[0][smallest]), (float)(sign * basis[1][smallest]),
                          (float)(sign * basis[2][smallest])};
}

typedef struct normal_edge {
    uint32_t from, to;
    float weight;
    bool flip;
} normal_edge;

typedef struct normal_set {
    uint32_t root;
    uint8_t rank;
    bool flip;
} normal_set;

static uint32_t normal_root(normal_set *sets, uint32_t vertex) {
    uint32_t parent = sets[vertex].root;
    if (parent != vertex) {
        sets[vertex].root = normal_root(sets, parent);
        sets[vertex].flip ^= sets[parent].flip;
    }
    return sets[vertex].root;
}

static int compare_normal_edges(const void *left, const void *right) {
    const normal_edge *a = left, *b = right;
    if (a->weight != b->weight) return a->weight < b->weight ? -1 : 1;
    if (a->from != b->from) return a->from < b->from ? -1 : 1;
    return (a->to > b->to) - (a->to < b->to);
}

int geomlib_pointcloud_estimate_normals(geomlib_pointcloud *pc, int k) {
    if (!pc || !pc->points || pc->point_count < 3 || pc->point_count > UINT32_MAX || k < 2 ||
        pc->point_capacity < pc->point_count) return GEOMLIB_MESH_ERROR_INVALID;
    size_t count = pc->point_count;
    size_t neighbors = (size_t)k < count ? (size_t)k : count - 1;
    if (pc->point_capacity > SIZE_MAX / sizeof(geomlib_vec3) ||
        count > SIZE_MAX / sizeof(normal_edge) / neighbors ||
        neighbors > SIZE_MAX / sizeof(normal_neighbor)) return GEOMLIB_MESH_ERROR_MEMORY;
    for (size_t i = 0; i < count; i++)
        if (!finite_point(pc->points[i])) return GEOMLIB_MESH_ERROR_INVALID;

    geomlib_vec3 *normals = calloc(pc->point_capacity, sizeof(geomlib_vec3));
    normal_neighbor *nearest = malloc(neighbors * sizeof(*nearest));
    normal_edge *edges = malloc(count * neighbors * sizeof(*edges));
    normal_set *sets = calloc(count, sizeof(*sets));
    int status = GEOMLIB_MESH_ERROR_MEMORY;
    if (!normals || !nearest || !edges || !sets) goto finish;

    size_t edge_count = 0;
    for (size_t i = 0; i < count; i++) {
        geomlib_vec3d point = point_as_double(pc->points[i]);
        size_t used = 0;
        for (size_t j = 0; j < count; j++) {
            if (j == i) continue;
            geomlib_vec3d delta = geomlib_vec3d_sub(point_as_double(pc->points[j]), point);
            double distance = geomlib_vec3d_dot(delta, delta);
            if (used == neighbors && distance >= nearest[used - 1].distance) continue;
            size_t at = used < neighbors ? used++ : used - 1;
            while (at > 0 && distance < nearest[at - 1].distance) {
                nearest[at] = nearest[at - 1];
                at--;
            }
            nearest[at] = (normal_neighbor){(uint32_t)j, distance};
        }
        normals[i] = estimate_pca_normal(pc, i, nearest, used);
        sets[i].root = (uint32_t)i;
        if (geomlib_vec3_dot(normals[i], normals[i]) == 0) continue;
        for (size_t j = 0; j < used; j++)
            edges[edge_count++] = (normal_edge){(uint32_t)i, nearest[j].index, 0, false};
    }

    size_t valid = 0;
    for (size_t i = 0; i < edge_count; i++) {
        normal_edge edge = edges[i];
        geomlib_vec3 target = normals[edge.to];
        if (geomlib_vec3_dot(target, target) == 0) continue;
        float alignment = geomlib_vec3_dot(normals[edge.from], target);
        edge.weight = 1 - fabsf(alignment);
        edge.flip = alignment < 0;
        edges[valid++] = edge;
    }
    qsort(edges, valid, sizeof(*edges), compare_normal_edges);
    /* Kruskal's forest aligns neighboring normals while preserving arbitrary global sign. */
    for (size_t i = 0; i < valid; i++) {
        normal_edge edge = edges[i];
        uint32_t a = normal_root(sets, edge.from), b = normal_root(sets, edge.to);
        if (a == b) continue;
        if (sets[a].rank < sets[b].rank) {
            uint32_t temp = a;
            a = b;
            b = temp;
        }
        sets[b].root = a;
        sets[b].flip = sets[edge.from].flip ^ sets[edge.to].flip ^ edge.flip;
        if (sets[a].rank == sets[b].rank) sets[a].rank++;
    }
    for (size_t i = 0; i < count; i++) {
        normal_root(sets, (uint32_t)i);
        if (sets[i].flip) normals[i] = geomlib_vec3_mul(normals[i], -1);
    }
    free(pc->normals);
    pc->normals = normals;
    normals = NULL;
    status = GEOMLIB_MESH_SUCCESS;

finish:
    free(normals);
    free(nearest);
    free(edges);
    free(sets);
    return status;
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
    (void)rotations;
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

/* ============================================================================
 * Shared reconstruction helpers
 * ============================================================================ */

typedef struct {
    uint32_t v[3];
} recon_face;

typedef struct {
    recon_face *faces;
    size_t count;
    size_t capacity;
} recon_face_list;

static void sort_triangle(uint32_t v[3]) {
    if (v[0] > v[1]) { uint32_t t = v[0]; v[0] = v[1]; v[1] = t; }
    if (v[1] > v[2]) { uint32_t t = v[1]; v[1] = v[2]; v[2] = t; }
    if (v[0] > v[1]) { uint32_t t = v[0]; v[0] = v[1]; v[1] = t; }
}

static int compare_recon_faces(const void *left, const void *right) {
    const recon_face *a = left, *b = right;
    uint32_t ka[3] = {a->v[0], a->v[1], a->v[2]};
    uint32_t kb[3] = {b->v[0], b->v[1], b->v[2]};
    sort_triangle(ka);
    sort_triangle(kb);
    for (int i = 0; i < 3; i++)
        if (ka[i] != kb[i]) return ka[i] < kb[i] ? -1 : 1;
    return 0;
}

static uint64_t edge_pair_key(uint32_t a, uint32_t b) {
    if (a > b) { uint32_t t = a; a = b; b = t; }
    return ((uint64_t)a << 32) | b;
}

typedef struct {
    uint64_t key;
    uint32_t face;
} recon_edge_record;

static int compare_recon_edges(const void *left, const void *right) {
    const recon_edge_record *a = left, *b = right;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return (a->face > b->face) - (a->face < b->face);
}

/* Deduplicate faces and keep every undirected edge at most twice. */
static void cleanup_faces(recon_face *faces, size_t *count) {
    size_t n = *count;
    if (n < 2) return;
    qsort(faces, n, sizeof(*faces), compare_recon_faces);
    size_t unique = 0;
    for (size_t i = 0; i < n; i++) {
        if (unique == 0 || compare_recon_faces(&faces[i], &faces[unique - 1]) != 0)
            faces[unique++] = faces[i];
    }
    n = unique;
    recon_edge_record *records = malloc(n * 3 * sizeof(*records));
    bool *dropped = calloc(n, sizeof(*dropped));
    if (!records || !dropped) {
        free(records);
        free(dropped);
        *count = n;
        return;
    }
    size_t entries = 0;
    for (size_t i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++) {
            records[entries].key = edge_pair_key(faces[i].v[j], faces[i].v[(j + 1) % 3]);
            records[entries].face = (uint32_t)i;
            entries++;
        }
    }
    qsort(records, entries, sizeof(*records), compare_recon_edges);
    for (size_t start = 0; start < entries;) {
        size_t end = start + 1;
        while (end < entries && records[end].key == records[start].key) end++;
        if (end - start > 2)
            for (size_t i = start + 2; i < end; i++) dropped[records[i].face] = true;
        start = end;
    }
    size_t write = 0;
    for (size_t i = 0; i < n; i++)
        if (!dropped[i]) faces[write++] = faces[i];
    *count = write;
    free(records);
    free(dropped);
}

static int face_list_push(recon_face_list *list, uint32_t a, uint32_t b, uint32_t c) {
    if (list->count == list->capacity) {
        size_t capacity = list->capacity ? list->capacity * 2 : 64;
        if (capacity > SIZE_MAX / sizeof(*list->faces)) return GEOMLIB_MESH_ERROR_MEMORY;
        recon_face *faces = realloc(list->faces, capacity * sizeof(*faces));
        if (!faces) return GEOMLIB_MESH_ERROR_MEMORY;
        list->faces = faces;
        list->capacity = capacity;
    }
    list->faces[list->count].v[0] = a;
    list->faces[list->count].v[1] = b;
    list->faces[list->count].v[2] = c;
    list->count++;
    return GEOMLIB_MESH_SUCCESS;
}

/* Uniform grid over the point cloud for radius queries. */
typedef struct recon_grid {
    size_t *starts;
    uint32_t *items;
    int64_t dim[3];
    geomlib_vec3 min;
    float cell;
} recon_grid;

static void recon_grid_free(recon_grid *grid) {
    if (!grid) return;
    free(grid->starts);
    free(grid->items);
    memset(grid, 0, sizeof(*grid));
}

static void recon_grid_index(const recon_grid *grid, geomlib_vec3 point, int64_t index[3]) {
    for (int i = 0; i < 3; i++) {
        float coordinate = ((&point.x)[i] - (&grid->min.x)[i]) / grid->cell;
        int64_t value = (int64_t)floorf(coordinate);
        if (value < 0) value = 0;
        if (value >= grid->dim[i]) value = grid->dim[i] - 1;
        index[i] = value;
    }
}

static int recon_grid_build(recon_grid *grid, const geomlib_pointcloud *pc, float cell) {
    memset(grid, 0, sizeof(*grid));
    if (!pc->point_count) return GEOMLIB_MESH_ERROR_INVALID;
    geomlib_aabb bounds = geomlib_pointcloud_bounds(pc);
    grid->min = bounds.min;
    grid->cell = fmaxf(cell, 1e-6f);
    uint64_t cells = 1;
    for (int i = 0; i < 3; i++) {
        float extent = (&bounds.max.x)[i] - (&bounds.min.x)[i];
        int64_t dim = (int64_t)floorf(fmaxf(extent, 0.0f) / grid->cell) + 1;
        if (dim < 1) dim = 1;
        if (dim > 1024) dim = 1024;
        grid->dim[i] = dim;
        if (cells > SIZE_MAX / (uint64_t)dim) return GEOMLIB_MESH_ERROR_MEMORY;
        cells *= (uint64_t)dim;
    }
    grid->starts = calloc((size_t)cells + 1, sizeof(*grid->starts));
    grid->items = malloc(pc->point_count * sizeof(*grid->items));
    if (!grid->starts || !grid->items) {
        recon_grid_free(grid);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t i = 0; i < pc->point_count; i++) {
        int64_t index[3];
        recon_grid_index(grid, pc->points[i], index);
        size_t cell_index = (size_t)((index[2] * grid->dim[1] + index[1]) * grid->dim[0] + index[0]);
        grid->starts[cell_index + 1]++;
    }
    for (size_t i = 0; i < (size_t)cells; i++) grid->starts[i + 1] += grid->starts[i];
    size_t *cursor = malloc((size_t)cells * sizeof(*cursor));
    if (!cursor) {
        recon_grid_free(grid);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    memcpy(cursor, grid->starts, (size_t)cells * sizeof(*cursor));
    for (size_t i = 0; i < pc->point_count; i++) {
        int64_t index[3];
        recon_grid_index(grid, pc->points[i], index);
        size_t cell_index = (size_t)((index[2] * grid->dim[1] + index[1]) * grid->dim[0] + index[0]);
        grid->items[cursor[cell_index]++] = (uint32_t)i;
    }
    free(cursor);
    return GEOMLIB_MESH_SUCCESS;
}

static size_t recon_grid_collect(const recon_grid *grid, const geomlib_pointcloud *pc,
                                 geomlib_vec3 center, float radius, uint32_t *out, size_t capacity) {
    int64_t lo[3], hi[3];
    for (int i = 0; i < 3; i++) {
        geomlib_vec3 offset = {0, 0, 0};
        (&offset.x)[i] = radius;
        int64_t low[3], high[3];
        recon_grid_index(grid, geomlib_vec3_sub(center, offset), low);
        recon_grid_index(grid, geomlib_vec3_add(center, offset), high);
        lo[i] = low[i];
        hi[i] = high[i];
    }
    size_t count = 0;
    float radius2 = radius * radius;
    for (int64_t z = lo[2]; z <= hi[2]; z++) {
        for (int64_t y = lo[1]; y <= hi[1]; y++) {
            for (int64_t x = lo[0]; x <= hi[0]; x++) {
                size_t cell_index = (size_t)((z * grid->dim[1] + y) * grid->dim[0] + x);
                for (size_t k = grid->starts[cell_index]; k < grid->starts[cell_index + 1]; k++) {
                    uint32_t index = grid->items[k];
                    geomlib_vec3 delta = geomlib_vec3_sub(pc->points[index], center);
                    if (geomlib_vec3_dot(delta, delta) > radius2) continue;
                    if (count < capacity) out[count] = index;
                    count++;
                }
            }
        }
    }
    return count;
}

static bool recon_grid_has_inside(const recon_grid *grid, const geomlib_pointcloud *pc,
                                  geomlib_vec3 center, float radius, const uint32_t exclude[3],
                                  size_t exclude_count) {
    int64_t lo[3], hi[3];
    for (int i = 0; i < 3; i++) {
        geomlib_vec3 offset = {0, 0, 0};
        (&offset.x)[i] = radius;
        int64_t low[3], high[3];
        recon_grid_index(grid, geomlib_vec3_sub(center, offset), low);
        recon_grid_index(grid, geomlib_vec3_add(center, offset), high);
        lo[i] = low[i];
        hi[i] = high[i];
    }
    float limit = radius * (1.0f - 1e-5f);
    float limit2 = limit * limit;
    for (int64_t z = lo[2]; z <= hi[2]; z++) {
        for (int64_t y = lo[1]; y <= hi[1]; y++) {
            for (int64_t x = lo[0]; x <= hi[0]; x++) {
                size_t cell_index = (size_t)((z * grid->dim[1] + y) * grid->dim[0] + x);
                for (size_t k = grid->starts[cell_index]; k < grid->starts[cell_index + 1]; k++) {
                    uint32_t index = grid->items[k];
                    bool skip = false;
                    for (size_t e = 0; e < exclude_count; e++) skip |= exclude[e] == index;
                    if (skip) continue;
                    geomlib_vec3 delta = geomlib_vec3_sub(pc->points[index], center);
                    if (geomlib_vec3_dot(delta, delta) < limit2) return true;
                }
            }
        }
    }
    return false;
}

static geomlib_vec3 estimate_point_normal(const geomlib_pointcloud *pc, size_t index,
                                          float radius) {
    geomlib_vec3 sum = {0, 0, 0};
    int pairs = 0;
    for (size_t j = 0; j < pc->point_count && pairs < 6; j++) {
        if (j == index) continue;
        geomlib_vec3 dj = geomlib_vec3_sub(pc->points[j], pc->points[index]);
        if (geomlib_vec3_length(dj) > radius) continue;
        for (size_t k = j + 1; k < pc->point_count && pairs < 6; k++) {
            if (k == index) continue;
            geomlib_vec3 dk = geomlib_vec3_sub(pc->points[k], pc->points[index]);
            if (geomlib_vec3_length(dk) > radius) continue;
            geomlib_vec3 n = geomlib_vec3_cross(dj, dk);
            if (geomlib_vec3_length(n) > 1e-12f) {
                sum = geomlib_vec3_add(sum, n);
                pairs++;
            }
        }
    }
    float length = geomlib_vec3_length(sum);
    return length > 1e-12f ? geomlib_vec3_mul(sum, 1.0f / length)
                           : geomlib_vec3_create(0, 0, 1);
}

/* Circumcenter/radius of a triangle in 3D (relative to p0). */
static geomlib_vec3 triangle_circumcenter_3d(geomlib_vec3 p0, geomlib_vec3 p1, geomlib_vec3 p2,
                                             float *radius) {
    geomlib_vec3d a = {p1.x - p0.x, p1.y - p0.y, p1.z - p0.z};
    geomlib_vec3d b = {p2.x - p0.x, p2.y - p0.y, p2.z - p0.z};
    double a2 = geomlib_vec3d_dot(a, a), b2 = geomlib_vec3d_dot(b, b), ab = geomlib_vec3d_dot(a, b);
    double det = a2 * b2 - ab * ab;
    if (!(fabs(det) > 1e-30)) {
        *radius = FLT_MAX;
        return p0;
    }
    double alpha = b2 * (a2 - ab) / (2.0 * det);
    double beta = a2 * (b2 - ab) / (2.0 * det);
    geomlib_vec3d offset = geomlib_vec3d_add(geomlib_vec3d_mul(a, alpha),
                                             geomlib_vec3d_mul(b, beta));
    *radius = (float)geomlib_vec3d_length(offset);
    return (geomlib_vec3){(float)(p0.x + offset.x), (float)(p0.y + offset.y),
                          (float)(p0.z + offset.z)};
}


static int reconstruct_greedy_projection(
    const geomlib_pointcloud *pc,
    const geomlib_surface_params *params,
    geomlib_trimesh *output
) {
    if (!pc || pc->point_count < 3) return GEOMLIB_MESH_ERROR_INVALID;
    float spacing = 1.0f;
    geomlib_pointcloud_get_stats(pc, &spacing, NULL);
    if (!(spacing > 0.0f)) spacing = 1.0f;
    float radius = params->search_radius > 0.0f ? params->search_radius : spacing * 3.0f;
    float max_edge = params->max_edge_length > 0.0f ? params->max_edge_length : spacing * 3.0f;
    float max_angle = params->max_surface_angle > 0.0f ? params->max_surface_angle : 45.0f;
    float cosine = cosf(max_angle * (float)M_PI / 180.0f);

    recon_grid grid;
    int status = recon_grid_build(&grid, pc, radius);
    if (status != GEOMLIB_MESH_SUCCESS) return status;

    recon_face_list faces = {0};
    uint32_t *ring = malloc(pc->point_count * sizeof(*ring));
    float *angles = malloc(pc->point_count * sizeof(*angles));
    if (!ring || !angles) {
        free(ring);
        free(angles);
        recon_grid_free(&grid);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t i = 0; i < pc->point_count && status == GEOMLIB_MESH_SUCCESS; i++) {
        geomlib_vec3 normal = pc->normals ? pc->normals[i]
                                          : estimate_point_normal(pc, i, radius);
        float length = geomlib_vec3_length(normal);
        normal = length > 1e-12f ? geomlib_vec3_mul(normal, 1.0f / length)
                                 : geomlib_vec3_create(0, 0, 1);
        geomlib_vec3 axis = fabsf(normal.x) < 0.9f ? geomlib_vec3_create(1, 0, 0)
                                                   : geomlib_vec3_create(0, 1, 0);
        geomlib_vec3 u = geomlib_vec3_normalize(
            geomlib_vec3_sub(axis, geomlib_vec3_mul(normal, geomlib_vec3_dot(axis, normal))));
        geomlib_vec3 w = geomlib_vec3_cross(normal, u);

        size_t count = recon_grid_collect(&grid, pc, pc->points[i], radius, ring,
                                          pc->point_count);
        if (count > pc->point_count) count = pc->point_count;
        size_t filtered = 0;
        for (size_t j = 0; j < count && filtered < pc->point_count; j++) {
            if (ring[j] == i) continue;
            geomlib_vec3 delta = geomlib_vec3_sub(pc->points[ring[j]], pc->points[i]);
            angles[filtered] = atan2f(geomlib_vec3_dot(delta, w), geomlib_vec3_dot(delta, u));
            ring[filtered++] = ring[j];
        }
        for (size_t a = 1; a < filtered; a++) {
            uint32_t key = ring[a];
            float angle = angles[a];
            size_t b = a;
            while (b > 0 && angles[b - 1] > angle) {
                angles[b] = angles[b - 1];
                ring[b] = ring[b - 1];
                b--;
            }
            angles[b] = angle;
            ring[b] = key;
        }
        for (size_t a = 0; a < filtered; a++) {
            uint32_t j = ring[a], k = ring[(a + 1) % filtered];
            if (j == k) continue;
            if (vec3_distance(pc->points[i], pc->points[j]) > max_edge ||
                vec3_distance(pc->points[i], pc->points[k]) > max_edge ||
                vec3_distance(pc->points[j], pc->points[k]) > max_edge) continue;
            geomlib_vec3 tri = geomlib_vec3_cross(
                geomlib_vec3_sub(pc->points[j], pc->points[i]),
                geomlib_vec3_sub(pc->points[k], pc->points[i]));
            float tri_length = geomlib_vec3_length(tri);
            if (!(tri_length > 1e-12f)) continue;
            tri = geomlib_vec3_mul(tri, 1.0f / tri_length);
            if (geomlib_vec3_dot(tri, normal) < cosine) continue;
            status = face_list_push(&faces, (uint32_t)i, j, k);
            if (status != GEOMLIB_MESH_SUCCESS) break;
        }
    }
    free(ring);
    free(angles);
    recon_grid_free(&grid);
    if (status != GEOMLIB_MESH_SUCCESS) {
        free(faces.faces);
        return status;
    }

    cleanup_faces(faces.faces, &faces.count);
    geomlib_trimesh_init(output, pc->point_count, faces.count);
    if ((pc->point_count && !output->vertices) || (faces.count && !output->indices)) {
        geomlib_trimesh_free(output);
        free(faces.faces);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t i = 0; i < pc->point_count; i++) {
        if (geomlib_trimesh_add_vertex(output, pc->points[i],
                pc->normals ? &pc->normals[i] : NULL, NULL) != GEOMLIB_MESH_SUCCESS) {
            geomlib_trimesh_free(output);
            free(faces.faces);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
    }
    for (size_t i = 0; i < faces.count; i++) {
        if (geomlib_trimesh_add_triangle(output, faces.faces[i].v[0], faces.faces[i].v[1],
                faces.faces[i].v[2]) != GEOMLIB_MESH_SUCCESS) {
            geomlib_trimesh_free(output);
            free(faces.faces);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
    }
    free(faces.faces);
    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Surface Reconstruction - Marching Cubes
 * ============================================================================ */

static int reconstruct_marching_cubes(
    const geomlib_pointcloud *pc,
    const geomlib_surface_params *params,
    geomlib_trimesh *output
) {
    if (!pc || pc->point_count < 3) return GEOMLIB_MESH_ERROR_INVALID;
    geomlib_aabb bounds = geomlib_pointcloud_bounds(pc);
    geomlib_vec3 extent = geomlib_vec3_sub(bounds.max, bounds.min);
    float spacing = 1.0f;
    geomlib_pointcloud_get_stats(pc, &spacing, NULL);
    if (!(spacing > 0.0f)) spacing = 1.0f;

    float voxel = params->voxel_size > 0.0f ? params->voxel_size : spacing;
    if (!(voxel > 0.0f) || !isfinite(voxel)) return GEOMLIB_MESH_ERROR_INVALID;
    float padding = 3.0f * voxel;
    for (;;) {
        double nx = ceil(((double)extent.x + 2.0 * padding) / voxel) + 1.0;
        double ny = ceil(((double)extent.y + 2.0 * padding) / voxel) + 1.0;
        double nz = ceil(((double)extent.z + 2.0 * padding) / voxel) + 1.0;
        if (nx <= 128.0 && ny <= 128.0 && nz <= 128.0) break;
        voxel *= 1.25f;
        padding = 3.0f * voxel;
        if (!isfinite(voxel)) return GEOMLIB_MESH_ERROR_INVALID;
    }
    uint32_t resolution[3] = {
        (uint32_t)(ceil(((double)extent.x + 2.0 * padding) / voxel) + 1.0),
        (uint32_t)(ceil(((double)extent.y + 2.0 * padding) / voxel) + 1.0),
        (uint32_t)(ceil(((double)extent.z + 2.0 * padding) / voxel) + 1.0)
    };
    if (resolution[0] < 2 || resolution[1] < 2 || resolution[2] < 2)
        return GEOMLIB_MESH_ERROR_INVALID;
    size_t total = (size_t)resolution[0] * resolution[1] * resolution[2];
    if (total > SIZE_MAX / sizeof(float)) return GEOMLIB_MESH_ERROR_MEMORY;
    float *values = malloc(total * sizeof(*values));
    if (!values) return GEOMLIB_MESH_ERROR_MEMORY;
    geomlib_vec3 origin = {bounds.min.x - padding, bounds.min.y - padding,
                           bounds.min.z - padding};

    /* Oriented-point implicit field: signed distance to the nearest sample plane. */
    recon_grid grid;
    int status = recon_grid_build(&grid, pc, fmaxf(2.0f * spacing, voxel));
    uint32_t *near = malloc(pc->point_count * sizeof(*near));
    if (status != GEOMLIB_MESH_SUCCESS || !near) {
        free(near);
        free(values);
        recon_grid_free(&grid);
        return status != GEOMLIB_MESH_SUCCESS ? status : GEOMLIB_MESH_ERROR_MEMORY;
    }
    float search = fmaxf(2.0f * spacing, voxel);
    float max_search = geomlib_vec3_length(extent) + 4.0f * padding;
    for (uint32_t z = 0; z < resolution[2]; z++) {
        for (uint32_t y = 0; y < resolution[1]; y++) {
            for (uint32_t x = 0; x < resolution[0]; x++) {
                geomlib_vec3 center = {
                    origin.x + ((float)x + 0.5f) * voxel,
                    origin.y + ((float)y + 0.5f) * voxel,
                    origin.z + ((float)z + 0.5f) * voxel
                };
                float radius = search;
                size_t count = 0;
                while (count == 0 && radius < max_search) {
                    count = recon_grid_collect(&grid, pc, center, radius, near,
                                               pc->point_count);
                    if (count > pc->point_count) count = pc->point_count;
                    if (!count) radius *= 2.0f;
                }
                float value = max_search;
                if (count) {
                    size_t nearest = near[0];
                    float best = vec3_distance(center, pc->points[nearest]);
                    for (size_t i = 1; i < count; i++) {
                        float distance = vec3_distance(center, pc->points[near[i]]);
                        if (distance < best) { best = distance; nearest = near[i]; }
                    }
                    geomlib_vec3 normal = pc->normals ? pc->normals[nearest]
                        : estimate_point_normal(pc, nearest, fmaxf(2.0f * spacing, voxel));
                    float length = geomlib_vec3_length(normal);
                    if (length > 1e-12f)
                        normal = geomlib_vec3_mul(normal, 1.0f / length);
                    value = geomlib_vec3_dot(geomlib_vec3_sub(center, pc->points[nearest]),
                                             normal);
                }
                values[((size_t)z * resolution[1] + y) * resolution[0] + x] = value;
            }
        }
    }
    free(near);
    recon_grid_free(&grid);

    /* Keep the level set away from the grid border so extraction is closed. */
    for (uint32_t z = 0; z < resolution[2]; z++) {
        for (uint32_t y = 0; y < resolution[1]; y++) {
            for (uint32_t x = 0; x < resolution[0]; x++) {
                if (x != 0 && y != 0 && z != 0 && x != resolution[0] - 1 &&
                    y != resolution[1] - 1 && z != resolution[2] - 1) continue;
                values[((size_t)z * resolution[1] + y) * resolution[0] + x] = max_search;
            }
        }
    }

    geomlib_sdf sdf;
    memset(&sdf, 0, sizeof(sdf));
    sdf.values = values;
    sdf.origin = origin;
    sdf.voxel_size = voxel;
    sdf.resolution[0] = resolution[0];
    sdf.resolution[1] = resolution[1];
    sdf.resolution[2] = resolution[2];
    sdf.total_voxels = total;
    sdf.extent = geomlib_vec3_create(resolution[0] * voxel, resolution[1] * voxel,
                                     resolution[2] * voxel);
    sdf.bounds.min = origin;
    sdf.bounds.max = geomlib_vec3_add(origin, sdf.extent);

    status = geomlib_sdf_to_mesh(&sdf, 0.0f, output);
    free(values);
    return status;
}

/* ============================================================================
 * Ball Pivoting Algorithm
 * ============================================================================ */

static int reconstruct_ball_pivoting(
    const geomlib_pointcloud *pc,
    const geomlib_surface_params *params,
    geomlib_trimesh *output
) {
    if (!pc || pc->point_count < 3) return GEOMLIB_MESH_ERROR_INVALID;
    float spacing = 1.0f;
    geomlib_pointcloud_get_stats(pc, &spacing, NULL);
    if (!(spacing > 0.0f)) spacing = 1.0f;
    float radius = params->ball_radius > 0.0f ? params->ball_radius : spacing * 2.0f;
    if (!(radius > 0.0f) || !isfinite(radius)) return GEOMLIB_MESH_ERROR_INVALID;

    recon_grid grid;
    int status = recon_grid_build(&grid, pc, 2.0f * radius);
    if (status != GEOMLIB_MESH_SUCCESS) return status;

    recon_face_list faces = {0};
    uint32_t *neighbors = malloc(pc->point_count * sizeof(*neighbors));
    if (!neighbors) {
        recon_grid_free(&grid);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    float radius2 = radius * radius;
    float pair_limit = 4.0f * radius2;
    for (size_t i = 0; i < pc->point_count && status == GEOMLIB_MESH_SUCCESS; i++) {
        size_t count = recon_grid_collect(&grid, pc, pc->points[i], 2.0f * radius, neighbors,
                                          pc->point_count);
        if (count > pc->point_count) count = pc->point_count;
        for (size_t a = 0; a < count && status == GEOMLIB_MESH_SUCCESS; a++) {
            uint32_t j = neighbors[a];
            if (j <= i) continue;
            geomlib_vec3 pj = pc->points[j];
            for (size_t b = a + 1; b < count; b++) {
                uint32_t k = neighbors[b];
                if (k <= j) continue;
                geomlib_vec3 pk = pc->points[k];
                if (geomlib_vec3_dot(geomlib_vec3_sub(pj, pk),
                                     geomlib_vec3_sub(pj, pk)) > pair_limit) continue;
                geomlib_vec3 triangle = geomlib_vec3_cross(
                    geomlib_vec3_sub(pj, pc->points[i]), geomlib_vec3_sub(pk, pc->points[i]));
                float area2 = geomlib_vec3_length(triangle);
                if (!(area2 > 1e-12f)) continue;
                geomlib_vec3 normal = geomlib_vec3_mul(triangle, 1.0f / area2);
                float circumradius;
                geomlib_vec3 center = triangle_circumcenter_3d(pc->points[i], pj, pk,
                                                               &circumradius);
                if (!(circumradius <= radius)) continue;
                bool flip = false;
                if (pc->normals) {
                    geomlib_vec3 average = geomlib_vec3_add(
                        geomlib_vec3_add(pc->normals[i], pc->normals[j]), pc->normals[k]);
                    float length = geomlib_vec3_length(average);
                    if (!(length > 1e-12f)) continue;
                    average = geomlib_vec3_mul(average, 1.0f / length);
                    float alignment = geomlib_vec3_dot(normal, average);
                    if (fabsf(alignment) < 0.2f) continue;
                    flip = alignment < 0.0f;
                }
                float height = sqrtf(fmaxf(radius2 - circumradius * circumradius, 0.0f));
                geomlib_vec3 offset = geomlib_vec3_mul(normal, height);
                uint32_t exclude[3] = {(uint32_t)i, j, k};
                bool empty = !recon_grid_has_inside(&grid, pc, geomlib_vec3_add(center, offset),
                                                     radius, exclude, 3) ||
                             !recon_grid_has_inside(&grid, pc, geomlib_vec3_sub(center, offset),
                                                     radius, exclude, 3);
                if (!empty) continue;
                status = flip ? face_list_push(&faces, (uint32_t)i, k, j)
                              : face_list_push(&faces, (uint32_t)i, j, k);
                if (status != GEOMLIB_MESH_SUCCESS) break;
            }
        }
    }
    free(neighbors);
    recon_grid_free(&grid);
    if (status != GEOMLIB_MESH_SUCCESS) {
        free(faces.faces);
        return status;
    }

    cleanup_faces(faces.faces, &faces.count);
    geomlib_trimesh_init(output, pc->point_count, faces.count);
    if ((pc->point_count && !output->vertices) || (faces.count && !output->indices)) {
        geomlib_trimesh_free(output);
        free(faces.faces);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t i = 0; i < pc->point_count; i++)
        if (geomlib_trimesh_add_vertex(output, pc->points[i],
                pc->normals ? &pc->normals[i] : NULL, NULL) != GEOMLIB_MESH_SUCCESS) {
            geomlib_trimesh_free(output);
            free(faces.faces);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
    for (size_t i = 0; i < faces.count; i++)
        if (geomlib_trimesh_add_triangle(output, faces.faces[i].v[0], faces.faces[i].v[1],
                faces.faces[i].v[2]) != GEOMLIB_MESH_SUCCESS) {
            geomlib_trimesh_free(output);
            free(faces.faces);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
    free(faces.faces);
    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Alpha Shape Algorithm
 * ============================================================================ */

typedef struct alpha_face_record {
    uint32_t key[3];
    recon_face face;
} alpha_face_record;

static int compare_alpha_faces(const void *left, const void *right) {
    const alpha_face_record *a = left, *b = right;
    for (int i = 0; i < 3; i++)
        if (a->key[i] != b->key[i]) return a->key[i] < b->key[i] ? -1 : 1;
    return 0;
}

static int reconstruct_alpha_shape(
    const geomlib_pointcloud *pc,
    const geomlib_surface_params *params,
    geomlib_trimesh *output
) {
    if (!pc || pc->point_count < 4) return GEOMLIB_MESH_ERROR_INVALID;
    geomlib_tetramesh tet = {0};
    int status = geomlib_tetrahedralize_points(pc->points, pc->point_count, &tet);
    if (status != GEOMLIB_MESH_SUCCESS) {
        if (status == GEOMLIB_MESH_ERROR_MEMORY) return status;
        return status == GEOMLIB_MESH_ERROR_DEGENERATE ? GEOMLIB_MESH_ERROR_DEGENERATE
                                                       : GEOMLIB_MESH_ERROR_INVALID;
    }
    float alpha = params->alpha;
    if (!(alpha > 0.0f)) {
        float spacing = 1.0f;
        geomlib_pointcloud_get_stats(pc, &spacing, NULL);
        alpha = spacing > 0.0f ? spacing * 1.5f : 1.0f;
    }
    size_t capacity = tet.tetrahedron_count * 4;
    alpha_face_record *records = malloc(capacity * sizeof(*records));
    if (!records) {
        geomlib_tetramesh_free(&tet);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    size_t count = 0;
    for (size_t cell = 0; cell < tet.tetrahedron_count; cell++) {
        const uint32_t *v = tet.tetrahedra[cell].vertices;
        geomlib_vec3 p[4] = {tet.vertices[v[0]], tet.vertices[v[1]], tet.vertices[v[2]],
                             tet.vertices[v[3]]};
        geomlib_vec3 center;
        float circumradius;
        geomlib_tetrahedron_circumsphere(p[0], p[1], p[2], p[3], &center, &circumradius);
        if (!(circumradius < alpha) || !isfinite(circumradius)) continue;
        for (int j = 0; j < 4; j++) {
            uint32_t face[3] = {v[(j + 1) % 4], v[(j + 2) % 4], v[(j + 3) % 4]};
            geomlib_vec3 a = tet.vertices[face[0]], b = tet.vertices[face[1]],
                          c = tet.vertices[face[2]];
            geomlib_vec3 normal = geomlib_vec3_cross(geomlib_vec3_sub(b, a),
                                                     geomlib_vec3_sub(c, a));
            if (geomlib_vec3_dot(normal, geomlib_vec3_sub(tet.vertices[v[j]], a)) > 0.0f) {
                uint32_t swap = face[1];
                face[1] = face[2];
                face[2] = swap;
            }
            records[count].face.v[0] = face[0];
            records[count].face.v[1] = face[1];
            records[count].face.v[2] = face[2];
            records[count].key[0] = face[0];
            records[count].key[1] = face[1];
            records[count].key[2] = face[2];
            sort_triangle(records[count].key);
            count++;
        }
    }
    qsort(records, count, sizeof(*records), compare_alpha_faces);
    recon_face_list faces = {0};
    for (size_t start = 0; start < count && status == GEOMLIB_MESH_SUCCESS;) {
        size_t end = start + 1;
        while (end < count && compare_alpha_faces(&records[start], &records[end]) == 0) end++;
        if (end - start == 1)
            status = face_list_push(&faces, records[start].face.v[0], records[start].face.v[1],
                                    records[start].face.v[2]);
        start = end;
    }
    free(records);
    geomlib_tetramesh_free(&tet);
    if (status != GEOMLIB_MESH_SUCCESS) {
        free(faces.faces);
        return status;
    }

    cleanup_faces(faces.faces, &faces.count);
    geomlib_trimesh_init(output, pc->point_count, faces.count);
    if ((pc->point_count && !output->vertices) || (faces.count && !output->indices)) {
        geomlib_trimesh_free(output);
        free(faces.faces);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t i = 0; i < pc->point_count; i++)
        if (geomlib_trimesh_add_vertex(output, pc->points[i],
                pc->normals ? &pc->normals[i] : NULL, NULL) != GEOMLIB_MESH_SUCCESS) {
            geomlib_trimesh_free(output);
            free(faces.faces);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
    for (size_t i = 0; i < faces.count; i++)
        if (geomlib_trimesh_add_triangle(output, faces.faces[i].v[0], faces.faces[i].v[1],
                faces.faces[i].v[2]) != GEOMLIB_MESH_SUCCESS) {
            geomlib_trimesh_free(output);
            free(faces.faces);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
    free(faces.faces);
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

    geomlib_surface_method method = target_faces
        ? GEOMLIB_SURFACE_MARCHING_CUBES
        : (pc->point_count < 1000 ? GEOMLIB_SURFACE_GREEDY_PROJECTION
                                  : GEOMLIB_SURFACE_MARCHING_CUBES);
    geomlib_surface_params params = geomlib_surface_params_default(method);
    if (target_faces) {
        float spacing = 0.0f;
        geomlib_pointcloud_get_stats(pc, &spacing, NULL);
        if (!(spacing > 0.0f)) spacing = 1.0f;
        double area = (double)pc->point_count * (double)spacing * (double)spacing;
        if (area > 0.0) {
            float edge = (float)sqrt(2.0 * area / (double)target_faces);
            float minimum = method == GEOMLIB_SURFACE_MARCHING_CUBES
                ? spacing * 0.5f : spacing * 1.5f;
            if (edge < minimum) edge = minimum;
            if (method == GEOMLIB_SURFACE_MARCHING_CUBES)
                params.voxel_size = edge;
            else {
                params.max_edge_length = edge;
                params.search_radius = edge * 2.0f;
            }
        }
    }
    return geomlib_reconstruct_surface(pc, &params, output);
}
