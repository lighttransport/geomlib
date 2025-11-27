/*
 * geomlib - Mesh Optimization Implementation
 */

#include "geomlib/geomlib_optimize.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

uint32_t geomlib_hash_position(geomlib_vec3 pos, float cell_size) {
    int32_t x = (int32_t)floorf(pos.x / cell_size);
    int32_t y = (int32_t)floorf(pos.y / cell_size);
    int32_t z = (int32_t)floorf(pos.z / cell_size);

    /* Simple hash function */
    uint32_t h = 2166136261u;
    h = (h ^ (uint32_t)x) * 16777619u;
    h = (h ^ (uint32_t)y) * 16777619u;
    h = (h ^ (uint32_t)z) * 16777619u;
    return h;
}

bool geomlib_vertices_equal(geomlib_vec3 a, geomlib_vec3 b, float tolerance) {
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    float dz = a.z - b.z;
    float dist_sq = dx*dx + dy*dy + dz*dz;
    return dist_sq <= tolerance * tolerance;
}

/* LRU cache for vertex cache simulation */
typedef struct vertex_cache {
    uint32_t *entries;
    uint32_t size;
    uint32_t count;
} vertex_cache;

static vertex_cache vertex_cache_create(uint32_t size) {
    vertex_cache cache;
    cache.entries = (uint32_t*)malloc(size * sizeof(uint32_t));
    cache.size = size;
    cache.count = 0;
    for (uint32_t i = 0; i < size; i++) {
        cache.entries[i] = UINT32_MAX;
    }
    return cache;
}

static void vertex_cache_free(vertex_cache *cache) {
    free(cache->entries);
    cache->entries = NULL;
}

static bool vertex_cache_lookup(vertex_cache *cache, uint32_t vertex) {
    /* Check if vertex is in cache */
    for (uint32_t i = 0; i < cache->count; i++) {
        if (cache->entries[i] == vertex) {
            /* Move to front (MRU) */
            if (i > 0) {
                uint32_t temp = cache->entries[i];
                memmove(&cache->entries[1], &cache->entries[0], i * sizeof(uint32_t));
                cache->entries[0] = temp;
            }
            return true;
        }
    }
    return false;
}

static void vertex_cache_insert(vertex_cache *cache, uint32_t vertex) {
    /* Insert at front, shift others back */
    if (cache->count < cache->size) {
        cache->count++;
    }
    memmove(&cache->entries[1], &cache->entries[0], (cache->count - 1) * sizeof(uint32_t));
    cache->entries[0] = vertex;
}

uint32_t geomlib_simulate_vertex_cache(
    const uint32_t *indices,
    size_t index_count,
    uint32_t cache_size
) {
    vertex_cache cache = vertex_cache_create(cache_size);
    uint32_t misses = 0;

    for (size_t i = 0; i < index_count; i++) {
        uint32_t vertex = indices[i];
        if (!vertex_cache_lookup(&cache, vertex)) {
            misses++;
            vertex_cache_insert(&cache, vertex);
        }
    }

    vertex_cache_free(&cache);
    return misses;
}

/* ============================================================================
 * Cache Analysis
 * ============================================================================ */

geomlib_cache_stats geomlib_analyze_vertex_cache(
    const geomlib_trimesh *mesh,
    uint32_t cache_size
) {
    geomlib_cache_stats stats;
    stats.cache_size = cache_size;

    /* Build index buffer */
    size_t index_count = mesh->triangle_count * 3;
    uint32_t *indices = (uint32_t*)malloc(index_count * sizeof(uint32_t));

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        indices[i * 3 + 0] = mesh->indices[i * 3 + 0];
        indices[i * 3 + 1] = mesh->indices[i * 3 + 1];
        indices[i * 3 + 2] = mesh->indices[i * 3 + 2];
    }

    /* Simulate cache */
    stats.cache_misses = geomlib_simulate_vertex_cache(indices, index_count, cache_size);
    stats.transformed_vertices = stats.cache_misses;

    /* Compute metrics */
    uint32_t triangle_count = (uint32_t)mesh->triangle_count;
    stats.acmr = (float)stats.cache_misses / (float)triangle_count;
    stats.atvr = (float)stats.transformed_vertices / (float)mesh->vertex_count;

    free(indices);
    return stats;
}

/* ============================================================================
 * Vertex Cache Optimization - Forsyth Algorithm
 * ============================================================================ */

#define FORSYTH_CACHE_SIZE_MAX 64
#define FORSYTH_VALENCE_BOOST_SCALE 2.0f
#define FORSYTH_VALENCE_BOOST_POWER 0.5f
#define FORSYTH_LAST_TRI_SCORE 0.75f
#define FORSYTH_CACHE_DECAY_POWER 1.5f

typedef struct vertex_data {
    uint32_t *triangles;        /* Triangles using this vertex */
    uint32_t triangle_count;
    uint32_t triangle_capacity;
    int32_t cache_position;     /* Position in cache (-1 if not in cache) */
    float score;
    bool emitted;
} vertex_data;

static float compute_vertex_score(const vertex_data *vertex, uint32_t cache_size) {
    if (vertex->triangle_count == 0) {
        return -1.0f;
    }

    float score = 0.0f;

    /* Cache position score */
    if (vertex->cache_position < 0) {
        /* Not in cache */
        score = 0.0f;
    } else if (vertex->cache_position < 3) {
        /* In last triangle - high score */
        score = FORSYTH_LAST_TRI_SCORE;
    } else {
        /* In cache but not in last triangle */
        float cache_decay = 1.0f - ((float)(vertex->cache_position - 3) / (float)(cache_size - 3));
        score = powf(cache_decay, FORSYTH_CACHE_DECAY_POWER);
    }

    /* Valence boost */
    float valence_boost = powf((float)vertex->triangle_count, -FORSYTH_VALENCE_BOOST_POWER);
    score += FORSYTH_VALENCE_BOOST_SCALE * valence_boost;

    return score;
}

static int geomlib_optimize_forsyth(
    const geomlib_trimesh *mesh,
    uint32_t *new_indices,
    uint32_t cache_size
) {
    size_t vertex_count = mesh->vertex_count;
    size_t triangle_count = mesh->triangle_count;

    /* Build vertex data */
    vertex_data *vertices = (vertex_data*)calloc(vertex_count, sizeof(vertex_data));

    /* Count triangles per vertex */
    for (size_t i = 0; i < triangle_count; i++) {
        for (int j = 0; j < 3; j++) {
            uint32_t v = mesh->indices[i * 3 + j];
            vertices[v].triangle_count++;
        }
    }

    /* Allocate triangle lists */
    for (size_t i = 0; i < vertex_count; i++) {
        vertices[i].triangle_capacity = vertices[i].triangle_count;
        vertices[i].triangles = (uint32_t*)malloc(vertices[i].triangle_capacity * sizeof(uint32_t));
        vertices[i].triangle_count = 0;
        vertices[i].cache_position = -1;
        vertices[i].emitted = false;
    }

    /* Build triangle lists */
    for (size_t i = 0; i < triangle_count; i++) {
        for (int j = 0; j < 3; j++) {
            uint32_t v = mesh->indices[i * 3 + j];
            vertices[v].triangles[vertices[v].triangle_count++] = (uint32_t)i;
        }
    }

    /* Initialize scores */
    for (size_t i = 0; i < vertex_count; i++) {
        vertices[i].score = compute_vertex_score(&vertices[i], cache_size);
    }

    /* Triangle scores and emitted flags */
    float *triangle_scores = (float*)malloc(triangle_count * sizeof(float));
    bool *triangle_emitted = (bool*)calloc(triangle_count, sizeof(bool));

    for (size_t i = 0; i < triangle_count; i++) {
        triangle_scores[i] = 0.0f;
        for (int j = 0; j < 3; j++) {
            uint32_t v = mesh->indices[i * 3 + j];
            triangle_scores[i] += vertices[v].score;
        }
    }

    /* LRU cache */
    uint32_t *cache = (uint32_t*)malloc(cache_size * sizeof(uint32_t));
    uint32_t cache_count = 0;
    for (uint32_t i = 0; i < cache_size; i++) {
        cache[i] = UINT32_MAX;
    }

    /* Emit triangles */
    size_t output_index = 0;
    size_t triangles_emitted = 0;

    while (triangles_emitted < triangle_count) {
        /* Find best triangle */
        uint32_t best_triangle = UINT32_MAX;
        float best_score = -1.0f;

        for (size_t i = 0; i < triangle_count; i++) {
            if (!triangle_emitted[i] && triangle_scores[i] > best_score) {
                best_score = triangle_scores[i];
                best_triangle = (uint32_t)i;
            }
        }

        if (best_triangle == UINT32_MAX) {
            /* No triangles with positive score, find any unemitted */
            for (size_t i = 0; i < triangle_count; i++) {
                if (!triangle_emitted[i]) {
                    best_triangle = (uint32_t)i;
                    break;
                }
            }
        }

        /* Emit triangle */
        triangle_emitted[best_triangle] = true;
        triangles_emitted++;

        uint32_t tri_vertices[3];
        for (int j = 0; j < 3; j++) {
            tri_vertices[j] = mesh->indices[best_triangle * 3 + j];
            new_indices[output_index++] = tri_vertices[j];
        }

        /* Update vertex cache and data */
        for (int j = 0; j < 3; j++) {
            uint32_t v = tri_vertices[j];

            /* Remove triangle from vertex list */
            for (uint32_t k = 0; k < vertices[v].triangle_count; k++) {
                if (vertices[v].triangles[k] == best_triangle) {
                    vertices[v].triangles[k] = vertices[v].triangles[vertices[v].triangle_count - 1];
                    vertices[v].triangle_count--;
                    break;
                }
            }

            /* Update cache position */
            int32_t cache_pos = -1;
            for (uint32_t k = 0; k < cache_count; k++) {
                if (cache[k] == v) {
                    cache_pos = (int32_t)k;
                    break;
                }
            }

            if (cache_pos >= 0) {
                /* Move to front */
                uint32_t temp = cache[cache_pos];
                memmove(&cache[1], &cache[0], cache_pos * sizeof(uint32_t));
                cache[0] = temp;
            } else {
                /* Insert at front */
                if (cache_count < cache_size) {
                    cache_count++;
                }
                memmove(&cache[1], &cache[0], (cache_count - 1) * sizeof(uint32_t));
                cache[0] = v;
            }
        }

        /* Update cache positions for all cached vertices */
        for (uint32_t i = 0; i < vertex_count; i++) {
            vertices[i].cache_position = -1;
        }
        for (uint32_t i = 0; i < cache_count; i++) {
            if (cache[i] != UINT32_MAX) {
                vertices[cache[i]].cache_position = (int32_t)i;
            }
        }

        /* Recompute scores for affected vertices and their triangles */
        for (int j = 0; j < 3; j++) {
            uint32_t v = tri_vertices[j];
            vertices[v].score = compute_vertex_score(&vertices[v], cache_size);

            /* Update triangle scores */
            for (uint32_t k = 0; k < vertices[v].triangle_count; k++) {
                uint32_t tri = vertices[v].triangles[k];
                if (!triangle_emitted[tri]) {
                    triangle_scores[tri] = 0.0f;
                    for (int m = 0; m < 3; m++) {
                        uint32_t tv = mesh->indices[tri * 3 + m];
                        triangle_scores[tri] += vertices[tv].score;
                    }
                }
            }
        }
    }

    /* Cleanup */
    for (size_t i = 0; i < vertex_count; i++) {
        free(vertices[i].triangles);
    }
    free(vertices);
    free(triangle_scores);
    free(triangle_emitted);
    free(cache);

    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_optimize_indices(
    const geomlib_trimesh *mesh,
    uint32_t *new_indices,
    geomlib_vcache_algorithm algorithm,
    uint32_t cache_size
) {
    if (!mesh || !new_indices) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    switch (algorithm) {
        case GEOMLIB_VCACHE_FORSYTH:
            return geomlib_optimize_forsyth(mesh, new_indices, cache_size);
        default:
            return GEOMLIB_MESH_ERROR_INVALID;
    }
}

int geomlib_optimize_vertex_cache(
    geomlib_trimesh *mesh,
    geomlib_vcache_algorithm algorithm,
    uint32_t cache_size
) {
    if (!mesh) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Optimize indices */
    size_t index_count = mesh->triangle_count * 3;
    uint32_t *new_indices = (uint32_t*)malloc(index_count * sizeof(uint32_t));

    int result = geomlib_optimize_indices(mesh, new_indices, algorithm, cache_size);
    if (result != GEOMLIB_MESH_SUCCESS) {
        free(new_indices);
        return result;
    }

    /* Update mesh faces */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        mesh->indices[i * 3 + 0] = new_indices[i * 3 + 0];
        mesh->indices[i * 3 + 1] = new_indices[i * 3 + 1];
        mesh->indices[i * 3 + 2] = new_indices[i * 3 + 2];
    }

    free(new_indices);
    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Spatial Hashing for Vertex Welding
 * ============================================================================ */

geomlib_spatial_hash geomlib_spatial_hash_create(
    const geomlib_vec3 *vertices,
    size_t vertex_count,
    float cell_size
) {
    geomlib_spatial_hash hash;

    /* Compute bounds */
    hash.min_bound = vertices[0];
    hash.max_bound = vertices[0];

    for (size_t i = 1; i < vertex_count; i++) {
        if (vertices[i].x < hash.min_bound.x) hash.min_bound.x = vertices[i].x;
        if (vertices[i].y < hash.min_bound.y) hash.min_bound.y = vertices[i].y;
        if (vertices[i].z < hash.min_bound.z) hash.min_bound.z = vertices[i].z;
        if (vertices[i].x > hash.max_bound.x) hash.max_bound.x = vertices[i].x;
        if (vertices[i].y > hash.max_bound.y) hash.max_bound.y = vertices[i].y;
        if (vertices[i].z > hash.max_bound.z) hash.max_bound.z = vertices[i].z;
    }

    hash.cell_size = cell_size;
    hash.grid_size = (uint32_t)(vertex_count / 4 + 1);  /* Reasonable hash table size */

    /* Allocate buckets */
    hash.grid = (uint32_t**)calloc(hash.grid_size, sizeof(uint32_t*));
    hash.bucket_sizes = (uint32_t*)calloc(hash.grid_size, sizeof(uint32_t));
    hash.bucket_capacities = (uint32_t*)calloc(hash.grid_size, sizeof(uint32_t));

    return hash;
}

void geomlib_spatial_hash_free(geomlib_spatial_hash *hash) {
    for (uint32_t i = 0; i < hash->grid_size; i++) {
        free(hash->grid[i]);
    }
    free(hash->grid);
    free(hash->bucket_sizes);
    free(hash->bucket_capacities);
}

static void spatial_hash_insert(geomlib_spatial_hash *hash, uint32_t vertex_idx, geomlib_vec3 pos) {
    uint32_t h = geomlib_hash_position(pos, hash->cell_size) % hash->grid_size;

    /* Ensure capacity */
    if (hash->bucket_sizes[h] >= hash->bucket_capacities[h]) {
        uint32_t new_capacity = hash->bucket_capacities[h] == 0 ? 4 : hash->bucket_capacities[h] * 2;
        hash->grid[h] = (uint32_t*)realloc(hash->grid[h], new_capacity * sizeof(uint32_t));
        hash->bucket_capacities[h] = new_capacity;
    }

    hash->grid[h][hash->bucket_sizes[h]++] = vertex_idx;
}

geomlib_weld_result geomlib_weld_vertices(
    const geomlib_trimesh *mesh,
    float tolerance
) {
    geomlib_weld_result result;
    result.original_vertex_count = mesh->vertex_count;

    /* Create spatial hash */
    geomlib_spatial_hash hash = geomlib_spatial_hash_create(
        mesh->vertices, mesh->vertex_count, tolerance * 2.0f
    );

    /* Build remap table */
    result.vertex_remap = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    bool *is_unique = (bool*)calloc(mesh->vertex_count, sizeof(bool));

    for (size_t i = 0; i < mesh->vertex_count; i++) {
        result.vertex_remap[i] = (uint32_t)i;
        is_unique[i] = true;
    }

    /* Insert vertices and find duplicates */
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        geomlib_vec3 pos = mesh->vertices[i];
        uint32_t h = geomlib_hash_position(pos, hash.cell_size) % hash.grid_size;

        /* Check bucket for nearby vertices */
        bool found_duplicate = false;
        for (uint32_t j = 0; j < hash.bucket_sizes[h]; j++) {
            uint32_t other_idx = hash.grid[h][j];
            if (other_idx < i && geomlib_vertices_equal(pos, mesh->vertices[other_idx], tolerance)) {
                /* Found duplicate */
                result.vertex_remap[i] = result.vertex_remap[other_idx];
                is_unique[i] = false;
                found_duplicate = true;
                break;
            }
        }

        if (!found_duplicate) {
            spatial_hash_insert(&hash, (uint32_t)i, pos);
        }
    }

    /* Count unique vertices */
    result.unique_vertex_count = 0;
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        if (is_unique[i]) {
            result.unique_vertex_count++;
        }
    }

    /* Build unique vertex list */
    result.vertices = (geomlib_vec3*)malloc(result.unique_vertex_count * sizeof(geomlib_vec3));
    uint32_t *compact_remap = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));

    uint32_t unique_idx = 0;
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        if (is_unique[i]) {
            result.vertices[unique_idx] = mesh->vertices[i];
            compact_remap[i] = unique_idx;
            unique_idx++;
        }
    }

    /* Update remap table to point to compact indices */
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        result.vertex_remap[i] = compact_remap[result.vertex_remap[i]];
    }

    /* Remap indices */
    size_t index_count = mesh->triangle_count * 3;
    result.indices = (uint32_t*)malloc(index_count * sizeof(uint32_t));

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        result.indices[i * 3 + 0] = result.vertex_remap[mesh->indices[i * 3 + 0]];
        result.indices[i * 3 + 1] = result.vertex_remap[mesh->indices[i * 3 + 1]];
        result.indices[i * 3 + 2] = result.vertex_remap[mesh->indices[i * 3 + 2]];
    }

    /* Cleanup */
    free(compact_remap);
    free(is_unique);
    geomlib_spatial_hash_free(&hash);

    return result;
}

int geomlib_apply_weld_result(
    geomlib_trimesh *mesh,
    const geomlib_weld_result *result
) {
    if (!mesh || !result) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Replace vertices */
    free(mesh->vertices);
    mesh->vertices = (geomlib_vec3*)malloc(result->unique_vertex_count * sizeof(geomlib_vec3));
    memcpy(mesh->vertices, result->vertices, result->unique_vertex_count * sizeof(geomlib_vec3));
    mesh->vertex_count = result->unique_vertex_count;

    /* Update indices */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        mesh->indices[i * 3 + 0] = result->indices[i * 3 + 0];
        mesh->indices[i * 3 + 1] = result->indices[i * 3 + 1];
        mesh->indices[i * 3 + 2] = result->indices[i * 3 + 2];
    }

    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_weld_result_free(geomlib_weld_result *result) {
    free(result->vertices);
    free(result->indices);
    free(result->vertex_remap);
}

/* ============================================================================
 * Vertex Fetch Optimization
 * ============================================================================ */

int geomlib_optimize_vertex_fetch(geomlib_trimesh *mesh) {
    if (!mesh) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Build vertex usage map */
    uint32_t *vertex_first_use = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    bool *vertex_used = (bool*)calloc(mesh->vertex_count, sizeof(bool));

    for (size_t i = 0; i < mesh->vertex_count; i++) {
        vertex_first_use[i] = UINT32_MAX;
    }

    /* Find first use of each vertex */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        for (int j = 0; j < 3; j++) {
            uint32_t v = mesh->indices[i * 3 + j];
            if (!vertex_used[v]) {
                vertex_first_use[v] = (uint32_t)i;
                vertex_used[v] = true;
            }
        }
    }

    /* Build new vertex order */
    uint32_t *new_order = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    uint32_t *remap = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));

    uint32_t new_idx = 0;
    for (uint32_t first_tri = 0; first_tri < mesh->triangle_count; first_tri++) {
        for (size_t v = 0; v < mesh->vertex_count; v++) {
            if (vertex_first_use[v] == first_tri) {
                new_order[new_idx] = (uint32_t)v;
                remap[v] = new_idx;
                new_idx++;
            }
        }
    }

    /* Reorder vertices */
    geomlib_vec3 *new_vertices = (geomlib_vec3*)malloc(mesh->vertex_count * sizeof(geomlib_vec3));
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        new_vertices[i] = mesh->vertices[new_order[i]];
    }

    free(mesh->vertices);
    mesh->vertices = new_vertices;

    /* Update indices */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        mesh->indices[i * 3 + 0] = remap[mesh->indices[i * 3 + 0]];
        mesh->indices[i * 3 + 1] = remap[mesh->indices[i * 3 + 1]];
        mesh->indices[i * 3 + 2] = remap[mesh->indices[i * 3 + 2]];
    }

    /* Cleanup */
    free(vertex_first_use);
    free(vertex_used);
    free(new_order);
    free(remap);

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Overdraw Optimization
 * ============================================================================ */

geomlib_overdraw_stats geomlib_analyze_overdraw(
    const geomlib_trimesh *mesh,
    geomlib_vec3 view_position,
    int resolution
) {
    geomlib_overdraw_stats stats;
    stats.avg_overdraw = 1.0f;  /* Placeholder */
    stats.max_overdraw = 1.0f;
    stats.total_shaded_pixels = 0;
    stats.visible_pixels = 0;
    return stats;
}

int geomlib_optimize_overdraw(
    geomlib_trimesh *mesh,
    geomlib_vec3 view_position
) {
    if (!mesh) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Compute triangle centers and depths */
    geomlib_vec3 *centers = (geomlib_vec3*)malloc(mesh->triangle_count * sizeof(geomlib_vec3));
    float *depths = (float*)malloc(mesh->triangle_count * sizeof(float));

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 v0 = mesh->vertices[mesh->indices[i * 3 + 0]];
        geomlib_vec3 v1 = mesh->vertices[mesh->indices[i * 3 + 1]];
        geomlib_vec3 v2 = mesh->vertices[mesh->indices[i * 3 + 2]];

        centers[i] = geomlib_vec3_create(
            (v0.x + v1.x + v2.x) / 3.0f,
            (v0.y + v1.y + v2.y) / 3.0f,
            (v0.z + v1.z + v2.z) / 3.0f
        );

        geomlib_vec3 to_view = geomlib_vec3_sub(view_position, centers[i]);
        depths[i] = geomlib_vec3_length(to_view);
    }

    /* Sort triangles front-to-back */
    uint32_t *sorted_indices = (uint32_t*)malloc(mesh->triangle_count * sizeof(uint32_t));
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        sorted_indices[i] = (uint32_t)i;
    }

    /* Simple bubble sort (could use qsort for larger meshes) */
    for (size_t i = 0; i < mesh->triangle_count - 1; i++) {
        for (size_t j = 0; j < mesh->triangle_count - i - 1; j++) {
            if (depths[sorted_indices[j]] > depths[sorted_indices[j + 1]]) {
                uint32_t temp = sorted_indices[j];
                sorted_indices[j] = sorted_indices[j + 1];
                sorted_indices[j + 1] = temp;
            }
        }
    }

    /* Reorder indices */
    uint32_t *new_indices = (uint32_t*)malloc(mesh->triangle_count * 3 * sizeof(uint32_t));
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t src_tri = sorted_indices[i];
        new_indices[i * 3 + 0] = mesh->indices[src_tri * 3 + 0];
        new_indices[i * 3 + 1] = mesh->indices[src_tri * 3 + 1];
        new_indices[i * 3 + 2] = mesh->indices[src_tri * 3 + 2];
    }

    free(mesh->indices);
    mesh->indices = new_indices;

    /* Cleanup */
    free(centers);
    free(depths);
    free(sorted_indices);

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Complete Optimization Pipeline
 * ============================================================================ */

geomlib_optimize_params geomlib_optimize_params_default(void) {
    geomlib_optimize_params params;
    params.flags = GEOMLIB_OPT_ALL;
    params.weld_tolerance = 1e-6f;
    params.cache_size = 16;
    params.vcache_algorithm = GEOMLIB_VCACHE_FORSYTH;
    params.view_position = geomlib_vec3_create(0, 0, 10);
    return params;
}

geomlib_optimize_result geomlib_optimize_mesh(
    geomlib_trimesh *mesh,
    const geomlib_optimize_params *params
) {
    geomlib_optimize_result result;
    result.vertices_before = mesh->vertex_count;

    clock_t start = clock();

    /* Analyze before */
    result.cache_before = geomlib_analyze_vertex_cache(mesh, params->cache_size);

    /* Weld vertices */
    if (params->flags & GEOMLIB_OPT_WELD_VERTICES) {
        geomlib_weld_result weld = geomlib_weld_vertices(mesh, params->weld_tolerance);
        geomlib_apply_weld_result(mesh, &weld);
        geomlib_weld_result_free(&weld);
    }

    /* Optimize vertex cache */
    if (params->flags & GEOMLIB_OPT_VERTEX_CACHE) {
        geomlib_optimize_vertex_cache(mesh, params->vcache_algorithm, params->cache_size);
    }

    /* Optimize overdraw */
    if (params->flags & GEOMLIB_OPT_OVERDRAW) {
        geomlib_optimize_overdraw(mesh, params->view_position);
    }

    /* Optimize vertex fetch */
    if (params->flags & GEOMLIB_OPT_VERTEX_FETCH) {
        geomlib_optimize_vertex_fetch(mesh);
    }

    /* Analyze after */
    result.cache_after = geomlib_analyze_vertex_cache(mesh, params->cache_size);
    result.vertices_after = mesh->vertex_count;

    clock_t end = clock();
    result.optimization_time_ms = (float)(end - start) / CLOCKS_PER_SEC * 1000.0f;

    return result;
}

/* ============================================================================
 * Triangle Strips (Stub Implementation)
 * ============================================================================ */

geomlib_tristrip_group geomlib_generate_tristrips(
    const geomlib_trimesh *mesh,
    bool use_restart
) {
    geomlib_tristrip_group group;
    group.strips = NULL;
    group.strip_count = 0;
    return group;
}

void geomlib_tristrip_group_free(geomlib_tristrip_group *group) {
    for (size_t i = 0; i < group->strip_count; i++) {
        free(group->strips[i].indices);
    }
    free(group->strips);
}

/* ============================================================================
 * Index Compression (Stub Implementation)
 * ============================================================================ */

geomlib_compressed_indices geomlib_compress_indices(
    const uint32_t *indices,
    size_t index_count,
    uint32_t vertex_count
) {
    geomlib_compressed_indices compressed;
    compressed.data = NULL;
    compressed.compressed_size = 0;
    compressed.original_count = index_count;
    compressed.vertex_count = vertex_count;
    return compressed;
}

int geomlib_decompress_indices(
    const geomlib_compressed_indices *compressed,
    uint32_t *indices
) {
    return GEOMLIB_MESH_ERROR_INVALID;
}

void geomlib_compressed_indices_free(geomlib_compressed_indices *compressed) {
    free(compressed->data);
}
