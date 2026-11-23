/*
 * geomlib - Mesh Optimization
 * Cache efficiency, vertex welding, and memory layout optimization
 */

#ifndef GEOMLIB_OPTIMIZE_H
#define GEOMLIB_OPTIMIZE_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Cache Efficiency Metrics
 * ============================================================================ */

/* Mesh cache statistics */
typedef struct geomlib_cache_stats {
    float acmr;                     /* Average Cache Miss Ratio */
    float atvr;                     /* Average Transform to Vertex Ratio */
    uint32_t transformed_vertices;  /* Number of vertex transformations */
    uint32_t cache_misses;          /* Number of cache misses */
    uint32_t cache_size;            /* Simulated cache size */
} geomlib_cache_stats;

/* Compute cache statistics for mesh */
geomlib_cache_stats geomlib_analyze_vertex_cache(
    const geomlib_trimesh *mesh,
    uint32_t cache_size  /* Typical: 16-32 for GPU vertex cache */
);

/* ============================================================================
 * Vertex Cache Optimization
 * ============================================================================ */

/* Vertex cache optimization algorithm */
typedef enum geomlib_vcache_algorithm {
    GEOMLIB_VCACHE_FORSYTH,         /* Tom Forsyth's linear-speed algorithm */
    GEOMLIB_VCACHE_TIPSIFY,         /* Tipsify algorithm */
    GEOMLIB_VCACHE_SIMPLE           /* Simple greedy algorithm */
} geomlib_vcache_algorithm;

/* Optimize mesh for vertex cache efficiency */
int geomlib_optimize_vertex_cache(
    geomlib_trimesh *mesh,
    geomlib_vcache_algorithm algorithm,
    uint32_t cache_size  /* Typical: 16-32 */
);

/* Reorder indices only (non-destructive) */
int geomlib_optimize_indices(
    const geomlib_trimesh *mesh,
    uint32_t *new_indices,  /* Output: reordered indices */
    geomlib_vcache_algorithm algorithm,
    uint32_t cache_size
);

/* ============================================================================
 * Vertex Welding and Deduplication
 * ============================================================================ */

/* Vertex welding result */
typedef struct geomlib_weld_result {
    geomlib_vec3 *vertices;         /* Unique vertices */
    uint32_t *indices;              /* Remapped indices */
    uint32_t *vertex_remap;         /* Old index -> new index mapping */
    size_t unique_vertex_count;
    size_t original_vertex_count;
} geomlib_weld_result;

/* Spatial hash grid for fast vertex lookup */
typedef struct geomlib_spatial_hash {
    uint32_t **grid;                /* Hash grid buckets */
    uint32_t *bucket_sizes;
    uint32_t *bucket_capacities;
    uint32_t grid_size;
    geomlib_vec3 min_bound;
    geomlib_vec3 max_bound;
    float cell_size;
} geomlib_spatial_hash;

/* Initialize spatial hash for vertex welding */
geomlib_spatial_hash geomlib_spatial_hash_create(
    const geomlib_vec3 *vertices,
    size_t vertex_count,
    float cell_size
);

/* Free spatial hash */
void geomlib_spatial_hash_free(geomlib_spatial_hash *hash);

/* Weld duplicate vertices using spatial hashing */
geomlib_weld_result geomlib_weld_vertices(
    const geomlib_trimesh *mesh,
    float tolerance  /* Distance threshold for considering vertices identical */
);

/* Apply welding result to mesh */
int geomlib_apply_weld_result(
    geomlib_trimesh *mesh,
    const geomlib_weld_result *result
);

/* Free weld result */
void geomlib_weld_result_free(geomlib_weld_result *result);

/* ============================================================================
 * Overdraw Optimization
 * ============================================================================ */

/* Overdraw statistics */
typedef struct geomlib_overdraw_stats {
    float avg_overdraw;             /* Average overdraw factor */
    float max_overdraw;             /* Maximum overdraw at any pixel */
    uint32_t total_shaded_pixels;   /* Total pixels shaded */
    uint32_t visible_pixels;        /* Visible pixels after depth test */
} geomlib_overdraw_stats;

/* Compute overdraw stats by rasterizing the mesh orthographically toward its
 * center. Back faces are culled; fragments are counted per pixel. */
geomlib_overdraw_stats geomlib_analyze_overdraw(
    const geomlib_trimesh *mesh,
    geomlib_vec3 view_position,
    int resolution  /* Resolution for overdraw buffer */
);

/* Optimize mesh to reduce overdraw */
int geomlib_optimize_overdraw(
    geomlib_trimesh *mesh,
    geomlib_vec3 view_position  /* Primary view direction */
);

/* ============================================================================
 * Vertex Fetch Optimization
 * ============================================================================ */

/* Reorder vertices for better memory locality */
int geomlib_optimize_vertex_fetch(
    geomlib_trimesh *mesh
);

/* Reorder vertex storage by first use, preserving triangle order and attributes.
 * Unreferenced vertices follow referenced vertices in their original order. */
int geomlib_reorder_vertices(
    geomlib_trimesh *mesh,
    const uint32_t *triangle_order  /* Optional permutation used for the first-use scan */
);

/* ============================================================================
 * Stripification
 * ============================================================================ */

/* Triangle strip */
typedef struct geomlib_tristrip {
    uint32_t *indices;
    size_t index_count;
    bool uses_restart;              /* Uses primitive restart */
    uint32_t restart_index;         /* Primitive restart index */
} geomlib_tristrip;

/* Triangle strip group */
typedef struct geomlib_tristrip_group {
    geomlib_tristrip *strips;
    size_t strip_count;
} geomlib_tristrip_group;

/* Generate triangle strips from mesh */
geomlib_tristrip_group geomlib_generate_tristrips(
    const geomlib_trimesh *mesh,
    bool use_restart  /* Use primitive restart for connecting strips */
);

/* Free triangle strips */
void geomlib_tristrip_group_free(geomlib_tristrip_group *group);

/* ============================================================================
 * Index Buffer Compression
 * ============================================================================ */

/* Compressed index buffer (delta encoding) */
typedef struct geomlib_compressed_indices {
    uint8_t *data;
    size_t compressed_size;
    size_t original_count;
    uint32_t vertex_count;
} geomlib_compressed_indices;

/* Compress indices using delta encoding */
geomlib_compressed_indices geomlib_compress_indices(
    const uint32_t *indices,
    size_t index_count,
    uint32_t vertex_count
);

/* Decompress indices */
int geomlib_decompress_indices(
    const geomlib_compressed_indices *compressed,
    uint32_t *indices  /* Output buffer (size: original_count) */
);

/* Free compressed indices */
void geomlib_compressed_indices_free(geomlib_compressed_indices *compressed);

/* ============================================================================
 * Complete Mesh Optimization Pipeline
 * ============================================================================ */

/* Optimization flags */
typedef enum geomlib_optimize_flags {
    GEOMLIB_OPT_WELD_VERTICES       = 1 << 0,  /* Remove duplicate vertices */
    GEOMLIB_OPT_VERTEX_CACHE        = 1 << 1,  /* Optimize for vertex cache */
    GEOMLIB_OPT_OVERDRAW            = 1 << 2,  /* Reduce overdraw */
    GEOMLIB_OPT_VERTEX_FETCH        = 1 << 3,  /* Optimize vertex fetch */
    GEOMLIB_OPT_ALL                 = 0xFF     /* All optimizations */
} geomlib_optimize_flags;

/* Optimization parameters */
typedef struct geomlib_optimize_params {
    uint32_t flags;                 /* Optimization flags */
    float weld_tolerance;           /* Vertex welding tolerance */
    uint32_t cache_size;            /* Vertex cache size (16-32) */
    geomlib_vcache_algorithm vcache_algorithm;
    geomlib_vec3 view_position;     /* For overdraw optimization */
} geomlib_optimize_params;

/* Optimization result */
typedef struct geomlib_optimize_result {
    geomlib_cache_stats cache_before;
    geomlib_cache_stats cache_after;
    size_t vertices_before;
    size_t vertices_after;
    float optimization_time_ms;
} geomlib_optimize_result;

/* Default optimization parameters */
geomlib_optimize_params geomlib_optimize_params_default(void);

/* Optimize mesh using specified parameters */
geomlib_optimize_result geomlib_optimize_mesh(
    geomlib_trimesh *mesh,
    const geomlib_optimize_params *params
);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/* Compute hash for position */
uint32_t geomlib_hash_position(
    geomlib_vec3 pos,
    float cell_size
);

/* Check if two vertices are within tolerance */
bool geomlib_vertices_equal(
    geomlib_vec3 a,
    geomlib_vec3 b,
    float tolerance
);

/* Simulate vertex cache and count misses */
uint32_t geomlib_simulate_vertex_cache(
    const uint32_t *indices,
    size_t index_count,
    uint32_t cache_size
);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_OPTIMIZE_H */
