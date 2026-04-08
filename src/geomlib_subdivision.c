/*
 * geomlib - Subdivision Surfaces Implementation
 * Loop subdivision for triangle meshes
 */

#include "geomlib/geomlib_allocator.h"
#include "geomlib/geomlib_c11.h"
#include "geomlib/geomlib_subdivision.h"
#include "geomlib/geomlib_healing.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Default subdivision parameters */
geomlib_subdivision_params geomlib_subdivision_params_default(
    geomlib_subdivision_scheme scheme,
    int levels
) {
    geomlib_subdivision_params params;
    params.scheme = scheme;
    params.levels = levels;
    params.preserve_boundaries = true;
    params.project_to_limit = false;
    return params;
}

/* Compute Loop subdivision weight (beta) for vertex with valence n */
float geomlib_loop_beta(int n) {
    if (n == 3) {
        return 3.0f / 16.0f;
    } else {
        /* Warren's improved weights */
        float angle = 2.0f * M_PI / n;
        float beta = (1.0f / n) * (5.0f / 8.0f - powf((3.0f / 8.0f + 0.25f * cosf(angle)), 2.0f));
        return beta;
    }
}

/* Edge structure for subdivision */
typedef struct edge_key {
    uint32_t v0, v1;  /* Sorted vertex indices */
} edge_key;

typedef struct edge_data {
    uint32_t new_vertex;  /* Index of newly created midpoint vertex */
    bool used;
} edge_data;

/* Create edge key (sorted) */
static edge_key make_edge_key(uint32_t a, uint32_t b) {
    edge_key key;
    if (a < b) {
        key.v0 = a;
        key.v1 = b;
    } else {
        key.v0 = b;
        key.v1 = a;
    }
    return key;
}

/* Simple hash map for edges */
#define EDGE_MAP_SIZE 8192

typedef struct edge_map_entry {
    edge_key key;
    uint32_t new_vertex;
    struct edge_map_entry *next;
} edge_map_entry;

typedef struct edge_map {
    edge_map_entry *buckets[EDGE_MAP_SIZE];
} edge_map;

static uint32_t edge_hash(edge_key key) {
    return (key.v0 * 73856093u + key.v1 * 19349663u) % EDGE_MAP_SIZE;
}

static void edge_map_init(edge_map *map) {
    geomlib_memset(map->buckets, 0, sizeof(map->buckets));
}

static void edge_map_put(edge_map *map, edge_key key, uint32_t new_vertex) {
    uint32_t hash = edge_hash(key);
    edge_map_entry *entry = geomlib_malloc(sizeof(edge_map_entry));
    entry->key = key;
    entry->new_vertex = new_vertex;
    entry->next = map->buckets[hash];
    map->buckets[hash] = entry;
}

static bool edge_map_get(edge_map *map, edge_key key, uint32_t *new_vertex) {
    uint32_t hash = edge_hash(key);
    edge_map_entry *entry = map->buckets[hash];
    while (entry) {
        if (entry->key.v0 == key.v0 && entry->key.v1 == key.v1) {
            *new_vertex = entry->new_vertex;
            return true;
        }
        entry = entry->next;
    }
    return false;
}

static void edge_map_free(edge_map *map) {
    for (int i = 0; i < EDGE_MAP_SIZE; i++) {
        edge_map_entry *entry = map->buckets[i];
        while (entry) {
            edge_map_entry *next = entry->next;
            geomlib_free(entry);
            entry = next;
        }
    }
}

/* Loop subdivision (one level) */
int geomlib_subdivide_loop(
    const geomlib_trimesh *input,
    geomlib_trimesh *output
) {
    if (!input || !output) return GEOMLIB_MESH_ERROR_INVALID;

    /* Build half-edge structure for neighbor queries */
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, input);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    /* Calculate output mesh size:
     * - Each original vertex remains
     * - Each edge gets a new vertex (midpoint)
     * - Each triangle splits into 4 triangles
     */
    size_t new_vertex_count = input->vertex_count + (hemesh.halfedge_count / 2);
    size_t new_triangle_count = input->triangle_count * 4;

    /* Initialize output mesh */
    geomlib_trimesh_init(output, new_vertex_count, new_triangle_count);

    /* Step 1: Copy original vertices with updated positions (odd vertices) */
    for (size_t v = 0; v < input->vertex_count; v++) {
        /* Find valence and neighbors */
        int valence = 0;
        geomlib_vec3 neighbor_sum = {0.0f, 0.0f, 0.0f};

        for (size_t he = 0; he < hemesh.halfedge_count; he++) {
            if (hemesh.halfedges[he].vertex == v) {
                uint32_t next = hemesh.halfedges[he].next;
                if (next != UINT32_MAX) {
                    uint32_t neighbor = hemesh.halfedges[next].vertex;
                    neighbor_sum = geomlib_vec3_add(neighbor_sum, input->vertices[neighbor]);
                    valence++;
                }
            }
        }

        geomlib_vec3 new_pos;
        if (valence > 0) {
            float beta = geomlib_loop_beta(valence);
            float alpha = 1.0f - valence * beta;

            /* new_pos = alpha * old_pos + beta * sum(neighbors) */
            new_pos = geomlib_vec3_add(
                geomlib_vec3_mul(input->vertices[v], alpha),
                geomlib_vec3_mul(neighbor_sum, beta)
            );
        } else {
            new_pos = input->vertices[v];
        }

        geomlib_trimesh_add_vertex(output, new_pos, NULL, NULL);
    }

    /* Step 2: Create edge vertices (even vertices) */
    edge_map emap;
    edge_map_init(&emap);

    for (size_t t = 0; t < input->triangle_count; t++) {
        uint32_t i0 = input->indices[t * 3 + 0];
        uint32_t i1 = input->indices[t * 3 + 1];
        uint32_t i2 = input->indices[t * 3 + 2];

        /* Process each edge */
        edge_key edges[3] = {
            make_edge_key(i0, i1),
            make_edge_key(i1, i2),
            make_edge_key(i2, i0)
        };

        for (int e = 0; e < 3; e++) {
            uint32_t new_v;
            if (edge_map_get(&emap, edges[e], &new_v)) {
                continue;  /* Edge already processed */
            }

            /* Compute edge midpoint with Loop weights */
            geomlib_vec3 v0 = input->vertices[edges[e].v0];
            geomlib_vec3 v1 = input->vertices[edges[e].v1];

            /* For boundary edges, use simple midpoint */
            /* For interior edges, use 3/8 weight on endpoints, 1/8 on opposite vertices */
            geomlib_vec3 edge_pos;

            /* Find opposite vertices in the two triangles sharing this edge */
            /* Simplified: use 3/8 weights (standard Loop) */
            edge_pos = geomlib_vec3_mul(geomlib_vec3_add(v0, v1), 3.0f / 8.0f);

            /* Find opposite vertices */
            uint32_t opp0 = UINT32_MAX, opp1 = UINT32_MAX;
            for (size_t tt = 0; tt < input->triangle_count; tt++) {
                uint32_t ti0 = input->indices[tt * 3 + 0];
                uint32_t ti1 = input->indices[tt * 3 + 1];
                uint32_t ti2 = input->indices[tt * 3 + 2];

                bool has_v0 = (ti0 == edges[e].v0 || ti1 == edges[e].v0 || ti2 == edges[e].v0);
                bool has_v1 = (ti0 == edges[e].v1 || ti1 == edges[e].v1 || ti2 == edges[e].v1);

                if (has_v0 && has_v1) {
                    /* This triangle contains the edge */
                    if (ti0 != edges[e].v0 && ti0 != edges[e].v1) {
                        if (opp0 == UINT32_MAX) opp0 = ti0; else opp1 = ti0;
                    }
                    if (ti1 != edges[e].v0 && ti1 != edges[e].v1) {
                        if (opp0 == UINT32_MAX) opp0 = ti1; else opp1 = ti1;
                    }
                    if (ti2 != edges[e].v0 && ti2 != edges[e].v1) {
                        if (opp0 == UINT32_MAX) opp0 = ti2; else opp1 = ti2;
                    }
                }
            }

            /* Add opposite vertex contributions */
            if (opp0 != UINT32_MAX) {
                edge_pos = geomlib_vec3_add(edge_pos,
                    geomlib_vec3_mul(input->vertices[opp0], 1.0f / 8.0f));
            }
            if (opp1 != UINT32_MAX) {
                edge_pos = geomlib_vec3_add(edge_pos,
                    geomlib_vec3_mul(input->vertices[opp1], 1.0f / 8.0f));
            }

            new_v = output->vertex_count;
            geomlib_trimesh_add_vertex(output, edge_pos, NULL, NULL);
            edge_map_put(&emap, edges[e], new_v);
        }
    }

    /* Step 3: Create new triangles */
    for (size_t t = 0; t < input->triangle_count; t++) {
        uint32_t v0 = input->indices[t * 3 + 0];
        uint32_t v1 = input->indices[t * 3 + 1];
        uint32_t v2 = input->indices[t * 3 + 2];

        /* Get edge midpoint vertices */
        uint32_t e01, e12, e20;
        edge_map_get(&emap, make_edge_key(v0, v1), &e01);
        edge_map_get(&emap, make_edge_key(v1, v2), &e12);
        edge_map_get(&emap, make_edge_key(v2, v0), &e20);

        /* Create 4 new triangles */
        geomlib_trimesh_add_triangle(output, v0, e01, e20);   /* Corner 0 */
        geomlib_trimesh_add_triangle(output, v1, e12, e01);   /* Corner 1 */
        geomlib_trimesh_add_triangle(output, v2, e20, e12);   /* Corner 2 */
        geomlib_trimesh_add_triangle(output, e01, e12, e20);  /* Center */
    }

    edge_map_free(&emap);
    geomlib_halfedge_free(&hemesh);

    return GEOMLIB_MESH_SUCCESS;
}

/* Simple midpoint subdivision */
int geomlib_subdivide_midpoint(
    const geomlib_trimesh *input,
    geomlib_trimesh *output
) {
    if (!input || !output) return GEOMLIB_MESH_ERROR_INVALID;

    size_t new_vertex_count = input->vertex_count + (input->triangle_count * 3);
    size_t new_triangle_count = input->triangle_count * 4;

    geomlib_trimesh_init(output, new_vertex_count, new_triangle_count);

    /* Copy original vertices */
    for (size_t v = 0; v < input->vertex_count; v++) {
        geomlib_trimesh_add_vertex(output, input->vertices[v], NULL, NULL);
    }

    /* Create edge map */
    edge_map emap;
    edge_map_init(&emap);

    /* Add midpoints and triangles */
    for (size_t t = 0; t < input->triangle_count; t++) {
        uint32_t v0 = input->indices[t * 3 + 0];
        uint32_t v1 = input->indices[t * 3 + 1];
        uint32_t v2 = input->indices[t * 3 + 2];

        /* Get or create edge midpoints */
        uint32_t e01, e12, e20;

        edge_key key01 = make_edge_key(v0, v1);
        if (!edge_map_get(&emap, key01, &e01)) {
            geomlib_vec3 mid = geomlib_vec3_mul(
                geomlib_vec3_add(input->vertices[v0], input->vertices[v1]), 0.5f);
            e01 = output->vertex_count;
            geomlib_trimesh_add_vertex(output, mid, NULL, NULL);
            edge_map_put(&emap, key01, e01);
        }

        edge_key key12 = make_edge_key(v1, v2);
        if (!edge_map_get(&emap, key12, &e12)) {
            geomlib_vec3 mid = geomlib_vec3_mul(
                geomlib_vec3_add(input->vertices[v1], input->vertices[v2]), 0.5f);
            e12 = output->vertex_count;
            geomlib_trimesh_add_vertex(output, mid, NULL, NULL);
            edge_map_put(&emap, key12, e12);
        }

        edge_key key20 = make_edge_key(v2, v0);
        if (!edge_map_get(&emap, key20, &e20)) {
            geomlib_vec3 mid = geomlib_vec3_mul(
                geomlib_vec3_add(input->vertices[v2], input->vertices[v0]), 0.5f);
            e20 = output->vertex_count;
            geomlib_trimesh_add_vertex(output, mid, NULL, NULL);
            edge_map_put(&emap, key20, e20);
        }

        /* Create 4 triangles */
        geomlib_trimesh_add_triangle(output, v0, e01, e20);
        geomlib_trimesh_add_triangle(output, v1, e12, e01);
        geomlib_trimesh_add_triangle(output, v2, e20, e12);
        geomlib_trimesh_add_triangle(output, e01, e12, e20);
    }

    edge_map_free(&emap);
    return GEOMLIB_MESH_SUCCESS;
}

/* Butterfly subdivision - simplified version */
int geomlib_subdivide_butterfly(
    const geomlib_trimesh *input,
    geomlib_trimesh *output
) {
    /* For simplicity, use midpoint subdivision */
    /* Full butterfly would require more complex stencil */
    return geomlib_subdivide_midpoint(input, output);
}

/* Subdivide mesh multiple levels */
int geomlib_subdivide_mesh(
    const geomlib_trimesh *input,
    const geomlib_subdivision_params *params,
    geomlib_trimesh *output
) {
    if (!input || !params || !output) return GEOMLIB_MESH_ERROR_INVALID;
    if (params->levels < 1) return GEOMLIB_MESH_ERROR_INVALID;

    /* Temporary mesh for iterative subdivision */
    geomlib_trimesh temp1, temp2;
    geomlib_trimesh *current = (geomlib_trimesh*)input;
    geomlib_trimesh *next = &temp1;

    int result = GEOMLIB_MESH_SUCCESS;

    for (int level = 0; level < params->levels; level++) {
        /* Choose subdivision method */
        switch (params->scheme) {
            case GEOMLIB_SUBDIV_LOOP:
                result = geomlib_subdivide_loop(current, next);
                break;
            case GEOMLIB_SUBDIV_BUTTERFLY:
                result = geomlib_subdivide_butterfly(current, next);
                break;
            case GEOMLIB_SUBDIV_MIDPOINT:
                result = geomlib_subdivide_midpoint(current, next);
                break;
            default:
                result = GEOMLIB_MESH_ERROR_INVALID;
        }

        if (result != GEOMLIB_MESH_SUCCESS) {
            if (level > 0) geomlib_trimesh_free(current);
            return result;
        }

        /* Swap for next iteration */
        if (level > 0) {
            geomlib_trimesh_free(current);
        }
        current = next;
        next = (next == &temp1) ? &temp2 : &temp1;
    }

    /* Copy result to output */
    *output = *current;

    return GEOMLIB_MESH_SUCCESS;
}

/* Compute Loop limit position (simplified) */
geomlib_vec3 geomlib_loop_limit_position(
    const geomlib_trimesh *mesh,
    uint32_t vertex_idx
) {
    /* For simplicity, return vertex position */
    /* Full implementation would compute true limit surface */
    if (vertex_idx < mesh->vertex_count) {
        return mesh->vertices[vertex_idx];
    }
    return (geomlib_vec3){0.0f, 0.0f, 0.0f};
}

/* Compute Loop limit normal (simplified) */
geomlib_vec3 geomlib_loop_limit_normal(
    const geomlib_trimesh *mesh,
    uint32_t vertex_idx
) {
    /* Return computed normal if available */
    if (mesh->normals && vertex_idx < mesh->vertex_count) {
        return mesh->normals[vertex_idx];
    }
    return (geomlib_vec3){0.0f, 1.0f, 0.0f};
}
