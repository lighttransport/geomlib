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

/* Tipsify: traverse the triangle adjacency graph, staying in cache as long
 * as possible. Implements the standard live-triangle/FIFO strategy. */
static int geomlib_optimize_tipsify(
    const geomlib_trimesh *mesh,
    uint32_t *new_indices,
    uint32_t cache_size
) {
    size_t V = mesh->vertex_count, T = mesh->triangle_count;
    size_t *vtri_count, *vtri_fill, *vtri;
    uint32_t *live, *stamp;
    uint8_t *emitted;
    uint32_t *fifo;
    size_t nfifo = 0;
    uint32_t time = cache_size + 1;
    size_t written = 0, i, k;

    vtri_count = (size_t *)calloc(V + 1, sizeof(size_t));
    live = (uint32_t *)calloc(V ? V : 1, sizeof(uint32_t));
    stamp = (uint32_t *)calloc(V ? V : 1, sizeof(uint32_t));
    emitted = (uint8_t *)calloc(T ? T : 1, sizeof(uint8_t));
    fifo = (uint32_t *)malloc(3 * T * sizeof(uint32_t) + 4);
    vtri = (size_t *)malloc(3 * T * sizeof(size_t) + 4);
    if (!vtri_count || !live || !stamp || !emitted || !fifo || !vtri) {
        free(vtri_count); free(live); free(stamp); free(emitted); free(fifo); free(vtri);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < T; ++i) {
        for (k = 0; k < 3; ++k) vtri_count[mesh->indices[i*3+k] + 1]++;
    }
    for (i = 0; i < V; ++i) vtri_count[i+1] += vtri_count[i];
    vtri_fill = (size_t *)calloc(V ? V : 1, sizeof(size_t));
    if (!vtri_fill) {
        free(vtri_count); free(live); free(stamp); free(emitted); free(fifo); free(vtri);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < T; ++i) {
        for (k = 0; k < 3; ++k) {
            uint32_t v = mesh->indices[i*3+k];
            vtri[vtri_count[v] + vtri_fill[v]++] = i;
        }
    }
    for (i = 0; i < V; ++i) live[i] = (uint32_t)(vtri_count[i+1] - vtri_count[i]);

    while (written < T) {
        size_t current = (size_t)-1;
        /* Pick a live vertex as a restart point if the FIFO is empty. */
        if (nfifo == 0) {
            size_t best_v = (size_t)-1;
            uint32_t best_live = 0;
            for (i = 0; i < V; ++i) {
                if (live[i] > best_live) { best_live = live[i]; best_v = i; }
            }
            if (best_v != (size_t)-1) {
                for (i = vtri_count[best_v]; i < vtri_count[best_v+1]; ++i) {
                    if (!emitted[vtri[i]]) { current = vtri[i]; break; }
                }
            }
            if (current == (size_t)-1) {
                for (i = 0; i < T; ++i) if (!emitted[i]) { current = i; break; }
            }
        } else {
            while (nfifo > 0) {
                uint32_t v = fifo[--nfifo];
                if (stamp[v] + cache_size > time) {
                    /* Vertex is still cached: prefer a neighbour using it. */
                    size_t best = (size_t)-1;
                    int best_score = -1;
                    for (i = vtri_count[v]; i < vtri_count[v+1]; ++i) {
                        size_t t = vtri[i];
                        int score = 0;
                        if (emitted[t]) continue;
                        for (k = 0; k < 3; ++k) {
                            uint32_t w = mesh->indices[t*3+k];
                            if (stamp[w] + cache_size > time) score++;
                        }
                        if (score > best_score) { best_score = score; best = t; }
                    }
                    if (best != (size_t)-1) { current = best; break; }
                } else if (live[v] > 0) {
                    for (i = vtri_count[v]; i < vtri_count[v+1]; ++i) {
                        if (!emitted[vtri[i]]) { current = vtri[i]; break; }
                    }
                    if (current != (size_t)-1) break;
                }
            }
            if (current == (size_t)-1) {
                for (i = 0; i < T; ++i) if (!emitted[i]) { current = i; break; }
            }
        }

        emitted[current] = 1;
        for (k = 0; k < 3; ++k) {
            uint32_t v = mesh->indices[current*3+k];
            new_indices[written*3+k] = v;
            if (live[v] > 0) live[v]--;
            stamp[v] = time++;
            if (live[v] > 0 && nfifo < 3 * T) fifo[nfifo++] = v;
        }
        written++;
    }

    free(vtri_count); free(vtri_fill); free(vtri);
    free(live); free(stamp); free(emitted); free(fifo);
    return GEOMLIB_MESH_SUCCESS;
}

/* Simple greedy: at each step pick the unemitted triangle with the most
 * corners currently in the cache. O(T^2) but deterministic and robust. */
static int geomlib_optimize_simple(
    const geomlib_trimesh *mesh,
    uint32_t *new_indices,
    uint32_t cache_size
) {
    size_t T = mesh->triangle_count, written = 0, i, k;
    uint8_t *emitted = (uint8_t *)calloc(T ? T : 1, sizeof(uint8_t));
    uint32_t *stamp = (uint32_t *)calloc(mesh->vertex_count ? mesh->vertex_count : 1,
                                         sizeof(uint32_t));
    uint32_t time = cache_size + 1;
    if (!emitted || !stamp) { free(emitted); free(stamp); return GEOMLIB_MESH_ERROR_MEMORY; }
    while (written < T) {
        size_t best = (size_t)-1;
        int best_score = -1;
        for (i = 0; i < T; ++i) {
            int score = 0;
            if (emitted[i]) continue;
            for (k = 0; k < 3; ++k) {
                uint32_t v = mesh->indices[i*3+k];
                if (stamp[v] + cache_size > time) score++;
            }
            if (score > best_score) { best_score = score; best = i; }
            if (score == 3) break;
        }
        if (best == (size_t)-1) break;
        emitted[best] = 1;
        for (k = 0; k < 3; ++k) {
            uint32_t v = mesh->indices[best*3+k];
            new_indices[written*3+k] = v;
            stamp[v] = time++;
        }
        written++;
    }
    free(emitted); free(stamp);
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
        case GEOMLIB_VCACHE_TIPSIFY:
            return geomlib_optimize_tipsify(mesh, new_indices, cache_size);
        case GEOMLIB_VCACHE_SIMPLE:
            return geomlib_optimize_simple(mesh, new_indices, cache_size);
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

int geomlib_reorder_vertices(geomlib_trimesh *mesh, const uint32_t *triangle_order) {
    if (!mesh || mesh->vertex_count > UINT32_MAX || mesh->triangle_count > UINT32_MAX ||
        (mesh->vertex_count && !mesh->vertices) ||
        (mesh->triangle_count && !mesh->indices) ||
        mesh->vertex_count > SIZE_MAX / sizeof(geomlib_vec3) ||
        mesh->triangle_count > SIZE_MAX / (3 * sizeof(uint32_t)))
        return GEOMLIB_MESH_ERROR_INVALID;
    if (!mesh->triangle_count) return GEOMLIB_MESH_SUCCESS;
    if (!mesh->vertex_count) return GEOMLIB_MESH_ERROR_INVALID;

    uint32_t *remap = malloc(mesh->vertex_count * sizeof(uint32_t));
    bool *seen = triangle_order ? calloc(mesh->triangle_count, sizeof(bool)) : NULL;
    geomlib_vec3 *vertices = malloc(mesh->vertex_count * sizeof(geomlib_vec3));
    geomlib_vec3 *normals = mesh->normals ? malloc(mesh->vertex_count * sizeof(geomlib_vec3)) : NULL;
    geomlib_vec2 *texcoords = mesh->texcoords ? malloc(mesh->vertex_count * sizeof(geomlib_vec2)) : NULL;
    int status = GEOMLIB_MESH_ERROR_MEMORY;
    if (!remap || (triangle_order && !seen) || !vertices ||
        (mesh->normals && !normals) || (mesh->texcoords && !texcoords)) goto finish;

    memset(remap, 0xff, mesh->vertex_count * sizeof(uint32_t));
    uint32_t next = 0;
    status = GEOMLIB_MESH_ERROR_INVALID;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        size_t triangle = triangle_order ? triangle_order[i] : i;
        if (triangle >= mesh->triangle_count || (seen && seen[triangle])) goto finish;
        if (seen) seen[triangle] = true;
        for (int j = 0; j < 3; j++) {
            uint32_t vertex = mesh->indices[triangle * 3 + j];
            if (vertex >= mesh->vertex_count) goto finish;
            if (remap[vertex] == UINT32_MAX) remap[vertex] = next++;
        }
    }

    for (size_t i = 0; i < mesh->vertex_count; i++) {
        if (remap[i] == UINT32_MAX) remap[i] = next++;
        vertices[remap[i]] = mesh->vertices[i];
        if (normals) normals[remap[i]] = mesh->normals[i];
        if (texcoords) texcoords[remap[i]] = mesh->texcoords[i];
    }
    for (size_t i = 0; i < mesh->triangle_count * 3; i++)
        mesh->indices[i] = remap[mesh->indices[i]];
    memcpy(mesh->vertices, vertices, mesh->vertex_count * sizeof(geomlib_vec3));
    if (normals) memcpy(mesh->normals, normals, mesh->vertex_count * sizeof(geomlib_vec3));
    if (texcoords) memcpy(mesh->texcoords, texcoords, mesh->vertex_count * sizeof(geomlib_vec2));
    status = GEOMLIB_MESH_SUCCESS;

finish:
    free(remap);
    free(seen);
    free(vertices);
    free(normals);
    free(texcoords);
    return status;
}

int geomlib_optimize_vertex_fetch(geomlib_trimesh *mesh) {
    return geomlib_reorder_vertices(mesh, NULL);
}

/* ============================================================================
 * Overdraw Optimization
 * ============================================================================ */

/* Orthographic projection toward the mesh center, with back-face culling and
 * per-pixel fragment counting. */
geomlib_overdraw_stats geomlib_analyze_overdraw(
    const geomlib_trimesh *mesh,
    geomlib_vec3 view_position,
    int resolution
) {
    geomlib_overdraw_stats stats = {0};
    if (!mesh || !mesh->vertices || !mesh->indices || !mesh->triangle_count ||
        !isfinite(view_position.x) || !isfinite(view_position.y) || !isfinite(view_position.z))
        return stats;
    if (resolution < 8) resolution = 8;
    if (resolution > 2048) resolution = 2048;

    geomlib_aabb bounds = geomlib_trimesh_compute_aabb(mesh);
    geomlib_vec3 center = geomlib_aabb_center(bounds);
    geomlib_vec3 forward = geomlib_vec3_sub(center, view_position);
    float forward_length = geomlib_vec3_length(forward);
    if (!(forward_length > 1e-9f)) forward = geomlib_vec3_create(0, 0, 1);
    forward = geomlib_vec3_normalize(forward);
    geomlib_vec3 up = fabsf(forward.y) < 0.9f ? geomlib_vec3_create(0, 1, 0)
                                             : geomlib_vec3_create(1, 0, 0);
    geomlib_vec3 right = geomlib_vec3_normalize(geomlib_vec3_cross(forward, up));
    up = geomlib_vec3_cross(right, forward);

    double min_x = DBL_MAX, max_x = -DBL_MAX, min_y = DBL_MAX, max_y = -DBL_MAX;
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        geomlib_vec3 delta = geomlib_vec3_sub(mesh->vertices[v], view_position);
        double x = geomlib_vec3_dot(delta, right);
        double y = geomlib_vec3_dot(delta, up);
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }
    double span = fmax(max_x - min_x, max_y - min_y);
    if (!(span > 1e-12)) return stats;
    double scale = (double)(resolution - 1) / span;
    double offset_x = (double)(resolution - 1) * 0.5 - 0.5 * (min_x + max_x) * scale;
    double offset_y = (double)(resolution - 1) * 0.5 - 0.5 * (min_y + max_y) * scale;

    size_t pixels = (size_t)resolution * (size_t)resolution;
    uint32_t *counts = calloc(pixels, sizeof(*counts));
    double *screen = malloc(mesh->vertex_count * 2 * sizeof(*screen));
    if (!counts || !screen) {
        free(counts);
        free(screen);
        return stats;
    }
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        geomlib_vec3 delta = geomlib_vec3_sub(mesh->vertices[v], view_position);
        screen[v * 2] = geomlib_vec3_dot(delta, right) * scale + offset_x;
        screen[v * 2 + 1] = geomlib_vec3_dot(delta, up) * scale + offset_y;
    }

    uint64_t fragments = 0;
    uint32_t max_count = 0;
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3], i1 = mesh->indices[t * 3 + 1],
                 i2 = mesh->indices[t * 3 + 2];
        if (i0 >= mesh->vertex_count || i1 >= mesh->vertex_count || i2 >= mesh->vertex_count)
            continue;
        geomlib_vec3 p0 = mesh->vertices[i0], p1 = mesh->vertices[i1], p2 = mesh->vertices[i2];
        geomlib_vec3 normal = geomlib_vec3_cross(geomlib_vec3_sub(p1, p0),
                                                 geomlib_vec3_sub(p2, p0));
        geomlib_vec3 centroid = geomlib_vec3_mul(
            geomlib_vec3_add(geomlib_vec3_add(p0, p1), p2), 1.0f / 3.0f);
        if (geomlib_vec3_dot(normal, geomlib_vec3_sub(view_position, centroid)) <= 0.0f)
            continue;
        double x0 = screen[i0 * 2], y0 = screen[i0 * 2 + 1];
        double x1 = screen[i1 * 2], y1 = screen[i1 * 2 + 1];
        double x2 = screen[i2 * 2], y2 = screen[i2 * 2 + 1];
        double area = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
        if (fabs(area) < 1e-12) continue;
        int lo_x = (int)floor(fmin(x0, fmin(x1, x2)));
        int hi_x = (int)ceil(fmax(x0, fmax(x1, x2)));
        int lo_y = (int)floor(fmin(y0, fmin(y1, y2)));
        int hi_y = (int)ceil(fmax(y0, fmax(y1, y2)));
        if (lo_x < 0) lo_x = 0;
        if (lo_y < 0) lo_y = 0;
        if (hi_x > resolution - 1) hi_x = resolution - 1;
        if (hi_y > resolution - 1) hi_y = resolution - 1;
        for (int y = lo_y; y <= hi_y; y++) {
            for (int x = lo_x; x <= hi_x; x++) {
                double px = x + 0.5, py = y + 0.5;
                double w0 = (x1 - px) * (y2 - py) - (y1 - py) * (x2 - px);
                double w1 = (x2 - px) * (y0 - py) - (y2 - py) * (x0 - px);
                double w2 = (x0 - px) * (y1 - py) - (y0 - py) * (x1 - px);
                if (area > 0.0) {
                    if (w0 < -1e-9 || w1 < -1e-9 || w2 < -1e-9) continue;
                } else {
                    if (w0 > 1e-9 || w1 > 1e-9 || w2 > 1e-9) continue;
                }
                size_t index = (size_t)y * resolution + (size_t)x;
                counts[index]++;
                fragments++;
                if (counts[index] > max_count) max_count = counts[index];
            }
        }
    }

    size_t visible = 0;
    for (size_t i = 0; i < pixels; i++) if (counts[i]) visible++;
    if (visible) stats.avg_overdraw = (float)((double)fragments / (double)visible);
    if (max_count) stats.max_overdraw = (float)max_count;
    stats.total_shaded_pixels = fragments > UINT32_MAX ? UINT32_MAX : (uint32_t)fragments;
    stats.visible_pixels = visible > UINT32_MAX ? UINT32_MAX : (uint32_t)visible;
    free(counts);
    free(screen);
    return stats;
}

typedef struct overdraw_depth_entry {
    uint32_t triangle;
    float depth;
} overdraw_depth_entry;

static int compare_overdraw_depth(const void *left, const void *right) {
    const overdraw_depth_entry *a = left, *b = right;
    return (a->depth > b->depth) - (a->depth < b->depth);
}

int geomlib_optimize_overdraw(
    geomlib_trimesh *mesh,
    geomlib_vec3 view_position
) {
    if (!mesh || !mesh->vertices || !mesh->indices || !mesh->triangle_count)
        return GEOMLIB_MESH_ERROR_INVALID;

    overdraw_depth_entry *entries = malloc(mesh->triangle_count * sizeof(*entries));
    if (!entries) return GEOMLIB_MESH_ERROR_MEMORY;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 v0 = mesh->vertices[mesh->indices[i * 3 + 0]];
        geomlib_vec3 v1 = mesh->vertices[mesh->indices[i * 3 + 1]];
        geomlib_vec3 v2 = mesh->vertices[mesh->indices[i * 3 + 2]];
        geomlib_vec3 center = geomlib_vec3_mul(
            geomlib_vec3_add(geomlib_vec3_add(v0, v1), v2), 1.0f / 3.0f);
        entries[i].triangle = (uint32_t)i;
        entries[i].depth = geomlib_vec3_length(geomlib_vec3_sub(view_position, center));
    }
    qsort(entries, mesh->triangle_count, sizeof(*entries), compare_overdraw_depth);

    uint32_t *new_indices = malloc(mesh->triangle_count * 3 * sizeof(*new_indices));
    if (!new_indices) {
        free(entries);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t source = entries[i].triangle;
        for (int j = 0; j < 3; j++)
            new_indices[i * 3 + j] = mesh->indices[source * 3 + j];
    }
    free(mesh->indices);
    mesh->indices = new_indices;
    free(entries);
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
 * Triangle Strips
 * ============================================================================ */

/* Build adjacency for stripification */
static void build_triangle_adjacency(
    const geomlib_trimesh *mesh,
    uint32_t **adjacency  /* 3 neighbors per triangle */
) {
    *adjacency = (uint32_t*)malloc(mesh->triangle_count * 3 * sizeof(uint32_t));

    /* Initialize all to no neighbor */
    for (size_t i = 0; i < mesh->triangle_count * 3; i++) {
        (*adjacency)[i] = UINT32_MAX;
    }

    /* Find neighbors by matching edges */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        for (int ei = 0; ei < 3; ei++) {
            if ((*adjacency)[i * 3 + ei] != UINT32_MAX) continue;

            uint32_t v0 = mesh->indices[i * 3 + ei];
            uint32_t v1 = mesh->indices[i * 3 + ((ei + 1) % 3)];

            /* Search for matching edge in other triangles */
            for (size_t j = i + 1; j < mesh->triangle_count; j++) {
                for (int ej = 0; ej < 3; ej++) {
                    uint32_t u0 = mesh->indices[j * 3 + ej];
                    uint32_t u1 = mesh->indices[j * 3 + ((ej + 1) % 3)];

                    /* Edges match if reversed (manifold mesh) */
                    if (v0 == u1 && v1 == u0) {
                        (*adjacency)[i * 3 + ei] = (uint32_t)j;
                        (*adjacency)[j * 3 + ej] = (uint32_t)i;
                        break;
                    }
                }
            }
        }
    }
}

/* Greedy strip extraction */
static geomlib_tristrip extract_strip(
    const geomlib_trimesh *mesh,
    const uint32_t *adjacency,
    bool *used,
    uint32_t start_tri,
    bool use_restart
) {
    geomlib_tristrip strip;
    strip.uses_restart = use_restart;
    strip.restart_index = UINT32_MAX;

    size_t capacity = 128;
    strip.indices = (uint32_t*)malloc(capacity * sizeof(uint32_t));
    strip.index_count = 0;

    /* Add first triangle */
    used[start_tri] = true;
    strip.indices[strip.index_count++] = mesh->indices[start_tri * 3 + 0];
    strip.indices[strip.index_count++] = mesh->indices[start_tri * 3 + 1];
    strip.indices[strip.index_count++] = mesh->indices[start_tri * 3 + 2];

    uint32_t current_tri = start_tri;

    /* Extend strip as far as possible */
    while (true) {
        /* Find next triangle */
        uint32_t next_tri = UINT32_MAX;
        uint32_t next_edge = 0;

        /* Try all edges */
        for (int e = 0; e < 3; e++) {
            uint32_t candidate = adjacency[current_tri * 3 + e];
            if (candidate != UINT32_MAX && !used[candidate]) {
                next_tri = candidate;
                next_edge = e;
                break;
            }
        }

        if (next_tri == UINT32_MAX) break;

        /* Add next vertex */
        used[next_tri] = true;

        /* Find the vertex not on the shared edge */
        uint32_t v0 = mesh->indices[current_tri * 3 + next_edge];
        uint32_t v1 = mesh->indices[current_tri * 3 + ((next_edge + 1) % 3)];

        for (int i = 0; i < 3; i++) {
            uint32_t v = mesh->indices[next_tri * 3 + i];
            if (v != v0 && v != v1) {
                if (strip.index_count >= capacity) {
                    capacity *= 2;
                    strip.indices = (uint32_t*)realloc(strip.indices, capacity * sizeof(uint32_t));
                }
                strip.indices[strip.index_count++] = v;
                break;
            }
        }

        current_tri = next_tri;
    }

    return strip;
}

geomlib_tristrip_group geomlib_generate_tristrips(
    const geomlib_trimesh *mesh,
    bool use_restart
) {
    geomlib_tristrip_group group;
    group.strips = NULL;
    group.strip_count = 0;

    if (!mesh || mesh->triangle_count == 0) {
        return group;
    }

    /* Build adjacency */
    uint32_t *adjacency;
    build_triangle_adjacency(mesh, &adjacency);

    /* Track used triangles */
    bool *used = (bool*)calloc(mesh->triangle_count, sizeof(bool));

    /* Extract strips */
    size_t capacity = 16;
    group.strips = (geomlib_tristrip*)malloc(capacity * sizeof(geomlib_tristrip));

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        if (used[i]) continue;

        geomlib_tristrip strip = extract_strip(mesh, adjacency, used, (uint32_t)i, use_restart);

        if (group.strip_count >= capacity) {
            capacity *= 2;
            group.strips = (geomlib_tristrip*)realloc(group.strips, capacity * sizeof(geomlib_tristrip));
        }

        group.strips[group.strip_count++] = strip;
    }

    free(adjacency);
    free(used);

    return group;
}

void geomlib_tristrip_group_free(geomlib_tristrip_group *group) {
    for (size_t i = 0; i < group->strip_count; i++) {
        free(group->strips[i].indices);
    }
    free(group->strips);
}

/* ============================================================================
 * Index Compression (Delta + Variable-Length Encoding)
 * ============================================================================ */

/* Encode unsigned int with variable-length encoding */
static size_t encode_varint(uint32_t value, uint8_t *buffer) {
    size_t bytes = 0;
    while (value >= 128) {
        buffer[bytes++] = (uint8_t)((value & 0x7F) | 0x80);
        value >>= 7;
    }
    buffer[bytes++] = (uint8_t)value;
    return bytes;
}

/* Decode variable-length encoded int */
static size_t decode_varint(const uint8_t *buffer, size_t length, uint32_t *value) {
    *value = 0;
    for (size_t byte_index = 0; byte_index < length && byte_index < 5; byte_index++) {
        uint32_t byte = buffer[byte_index];
        if (byte_index == 4 && byte > 15) return 0;
        *value |= (byte & 0x7f) << (7 * byte_index);
        if (byte < 128) return byte_index + 1;
    }
    return 0;
}

/* The delta is a signed 32-bit value represented modulo 2^32. */
static uint32_t zigzag_encode(uint32_t delta) {
    return (delta << 1) ^ (0u - (delta >> 31));
}

static uint32_t zigzag_decode(uint32_t value) {
    return (value >> 1) ^ (0u - (value & 1u));
}

geomlib_compressed_indices geomlib_compress_indices(
    const uint32_t *indices,
    size_t index_count,
    uint32_t vertex_count
) {
    geomlib_compressed_indices compressed = {0};
    compressed.original_count = index_count;
    compressed.vertex_count = vertex_count;

    if (!indices || index_count == 0 || index_count > SIZE_MAX / 5) {
        return compressed;
    }

    /* Worst case: five bytes per index. */
    size_t max_size = index_count * 5;
    compressed.data = (uint8_t*)malloc(max_size);
    if (!compressed.data) return compressed;

    /* Encode using delta encoding */
    uint32_t last_index = 0;

    for (size_t i = 0; i < index_count; i++) {
        /* Compute delta */
        uint32_t delta = indices[i] - last_index;

        /* Zig-zag encode delta (converts signed to unsigned) */
        uint32_t encoded = zigzag_encode(delta);

        /* Variable-length encode */
        compressed.compressed_size += encode_varint(
            encoded,
            compressed.data + compressed.compressed_size
        );

        last_index = indices[i];
    }

    /* Shrink to actual size */
    uint8_t *smaller = (uint8_t*)realloc(compressed.data, compressed.compressed_size);
    if (smaller) compressed.data = smaller;

    return compressed;
}

int geomlib_decompress_indices(
    const geomlib_compressed_indices *compressed,
    uint32_t *indices
) {
    if (!compressed || !compressed->data || !indices) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    size_t offset = 0;
    uint32_t last_index = 0;

    for (size_t i = 0; i < compressed->original_count; i++) {
        if (offset >= compressed->compressed_size) {
            return GEOMLIB_MESH_ERROR_INVALID;
        }

        /* Decode variable-length int */
        uint32_t encoded;
        size_t bytes = decode_varint(compressed->data + offset,
                                     compressed->compressed_size - offset, &encoded);
        if (!bytes) return GEOMLIB_MESH_ERROR_INVALID;
        offset += bytes;

        uint32_t delta = zigzag_decode(encoded);

        /* Reconstruct index */
        indices[i] = last_index + delta;
        last_index = indices[i];
    }

    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_compressed_indices_free(geomlib_compressed_indices *compressed) {
    free(compressed->data);
}
