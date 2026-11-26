/*
 * geomlib - Subdivision Surfaces Implementation
 * Loop subdivision for triangle meshes
 */

#include "geomlib/geomlib_subdivision.h"
#include "geomlib/geomlib_healing.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/* ============================================================================
 * Shared topology helpers for Butterfly and Loop limit evaluation
 * ============================================================================ */

typedef struct bf_raw_edge {
    uint32_t a, b;
    uint32_t face;
    uint32_t opposite;
} bf_raw_edge;

typedef struct bf_edge {
    uint32_t a, b;
    uint32_t face[2];
    uint32_t opposite[2];
    uint32_t new_vertex;
} bf_edge;

static int compare_raw_edges(const void *left, const void *right) {
    const bf_raw_edge *x = left, *y = right;
    if (x->a != y->a) return x->a < y->a ? -1 : 1;
    if (x->b != y->b) return x->b < y->b ? -1 : 1;
    return (x->face > y->face) - (x->face < y->face);
}

static void edge_pair(uint32_t x, uint32_t y, uint32_t *a, uint32_t *b) {
    if (x < y) { *a = x; *b = y; } else { *a = y; *b = x; }
}

static int edge_compare_key(uint32_t a0, uint32_t b0, uint32_t a1, uint32_t b1) {
    if (a0 != a1) return a0 < a1 ? -1 : 1;
    if (b0 != b1) return b0 < b1 ? -1 : 1;
    return 0;
}

/* Merge duplicated triangle edges into edges with up to two incident faces. */
static bf_edge *build_bf_edges(const geomlib_trimesh *mesh, size_t *edge_count) {
    if (mesh->triangle_count > SIZE_MAX / (3 * sizeof(bf_raw_edge))) return NULL;
    bf_raw_edge *raw = malloc(mesh->triangle_count * 3 * sizeof(*raw));
    bf_edge *edges = malloc(mesh->triangle_count * 3 * sizeof(*edges));
    if (!raw || !edges) {
        free(raw);
        free(edges);
        return NULL;
    }
    size_t count = 0;
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t v[3] = {mesh->indices[t * 3], mesh->indices[t * 3 + 1], mesh->indices[t * 3 + 2]};
        for (int j = 0; j < 3; j++) {
            bf_raw_edge entry;
            edge_pair(v[j], v[(j + 1) % 3], &entry.a, &entry.b);
            entry.face = (uint32_t)t;
            entry.opposite = v[(j + 2) % 3];
            raw[count++] = entry;
        }
    }
    qsort(raw, count, sizeof(*raw), compare_raw_edges);
    size_t unique = 0;
    for (size_t start = 0; start < count;) {
        size_t end = start + 1;
        while (end < count && raw[end].a == raw[start].a && raw[end].b == raw[start].b) end++;
        bf_edge edge;
        edge.a = raw[start].a;
        edge.b = raw[start].b;
        edge.face[0] = raw[start].face;
        edge.opposite[0] = raw[start].opposite;
        edge.face[1] = UINT32_MAX;
        edge.opposite[1] = UINT32_MAX;
        for (size_t i = start + 1; i < end; i++) {
            if (raw[i].face == edge.face[0]) continue;
            edge.face[1] = raw[i].face;
            edge.opposite[1] = raw[i].opposite;
            break;
        }
        edge.new_vertex = UINT32_MAX;
        edges[unique++] = edge;
        start = end;
    }
    free(raw);
    *edge_count = unique;
    return edges;
}

static bf_edge *find_bf_edge(bf_edge *edges, size_t count, uint32_t x, uint32_t y) {
    uint32_t a, b;
    edge_pair(x, y, &a, &b);
    size_t lo = 0, hi = count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int cmp = edge_compare_key(edges[mid].a, edges[mid].b, a, b);
        if (cmp == 0) return &edges[mid];
        if (cmp < 0) lo = mid + 1;
        else hi = mid;
    }
    return NULL;
}

/* Area-weighted vertex normals (unnormalized input faces). */
static geomlib_vec3 *area_normals(const geomlib_trimesh *mesh) {
    geomlib_vec3 *normals = calloc(mesh->vertex_count, sizeof(*normals));
    if (!normals) return NULL;
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t i0 = mesh->indices[t * 3], i1 = mesh->indices[t * 3 + 1],
                 i2 = mesh->indices[t * 3 + 2];
        geomlib_vec3 n = geomlib_vec3_cross(geomlib_vec3_sub(mesh->vertices[i1], mesh->vertices[i0]),
                                            geomlib_vec3_sub(mesh->vertices[i2], mesh->vertices[i0]));
        normals[i0] = geomlib_vec3_add(normals[i0], n);
        normals[i1] = geomlib_vec3_add(normals[i1], n);
        normals[i2] = geomlib_vec3_add(normals[i2], n);
    }
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        float length = geomlib_vec3_length(normals[v]);
        normals[v] = length > 1e-12f ? geomlib_vec3_mul(normals[v], 1.0f / length)
                                     : geomlib_vec3_create(0, 0, 1);
    }
    return normals;
}

static int compare_vertex_ids(const void *left, const void *right) {
    uint32_t a = *(const uint32_t*)left, b = *(const uint32_t*)right;
    return (a > b) - (a < b);
}

/* Unique one-ring neighbors sorted cyclically around each vertex normal. */
static int build_ordered_rings(const geomlib_trimesh *mesh, const geomlib_vec3 *normals,
                               size_t **offsets_out, uint32_t **neighbors_out,
                               size_t **counts_out) {
    size_t vertices = mesh->vertex_count;
    if (vertices > SIZE_MAX / (6 * sizeof(uint32_t)) ||
        mesh->triangle_count > SIZE_MAX / (6 * sizeof(uint32_t)))
        return GEOMLIB_MESH_ERROR_MEMORY;
    size_t *counts = calloc(vertices, sizeof(*counts));
    size_t *offsets = calloc(vertices + 1, sizeof(*offsets));
    uint32_t *neighbors = malloc(mesh->triangle_count * 6 * sizeof(*neighbors));
    if (!counts || !offsets || !neighbors) {
        free(counts);
        free(offsets);
        free(neighbors);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t t = 0; t < mesh->triangle_count; t++)
        for (int j = 0; j < 3; j++) {
            uint32_t v = mesh->indices[t * 3 + j];
            if (v < vertices) counts[v] += 2;
        }
    for (size_t i = 0; i < vertices; i++) offsets[i + 1] = offsets[i] + counts[i];
    size_t *cursor = malloc(vertices * sizeof(*cursor));
    if (!cursor) {
        free(counts); free(offsets); free(neighbors);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    memcpy(cursor, offsets, vertices * sizeof(*cursor));
    for (size_t t = 0; t < mesh->triangle_count; t++) {
        uint32_t v[3] = {mesh->indices[t * 3], mesh->indices[t * 3 + 1], mesh->indices[t * 3 + 2]};
        for (int j = 0; j < 3; j++) {
            if (v[j] >= vertices) continue;
            neighbors[cursor[v[j]]++] = v[(j + 1) % 3];
            neighbors[cursor[v[j]]++] = v[(j + 2) % 3];
        }
    }
    free(cursor);

    float *angles = malloc(offsets[vertices] * sizeof(*angles));
    if (!angles) {
        free(counts); free(offsets); free(neighbors);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t v = 0; v < vertices; v++) {
        size_t count = counts[v];
        uint32_t *start = neighbors + offsets[v];
        qsort(start, count, sizeof(*start), compare_vertex_ids);
        size_t unique = 0;
        for (size_t i = 0; i < count; i++)
            if (!unique || start[unique - 1] != start[i]) start[unique++] = start[i];
        counts[v] = unique;

        geomlib_vec3 n = normals[v];
        geomlib_vec3 axis = fabsf(n.x) < 0.9f ? geomlib_vec3_create(1, 0, 0)
                                              : geomlib_vec3_create(0, 1, 0);
        geomlib_vec3 u = geomlib_vec3_normalize(
            geomlib_vec3_sub(axis, geomlib_vec3_mul(n, geomlib_vec3_dot(axis, n))));
        geomlib_vec3 w = geomlib_vec3_cross(n, u);
        for (size_t i = 0; i < unique; i++) {
            geomlib_vec3 d = geomlib_vec3_sub(mesh->vertices[start[i]], mesh->vertices[v]);
            angles[i] = atan2f(geomlib_vec3_dot(d, w), geomlib_vec3_dot(d, u));
        }
        for (size_t i = 1; i < unique; i++) {
            uint32_t key = start[i];
            float angle = angles[i];
            size_t j = i;
            while (j > 0 && angles[j - 1] > angle) {
                angles[j] = angles[j - 1];
                start[j] = start[j - 1];
                j--;
            }
            angles[j] = angle;
            start[j] = key;
        }
    }
    free(angles);
    *offsets_out = offsets;
    *neighbors_out = neighbors;
    *counts_out = counts;
    return GEOMLIB_MESH_SUCCESS;
}

static bool edge_other_opposite(bf_edge *edges, size_t count, uint32_t x, uint32_t y,
                                uint32_t face, uint32_t *opposite) {
    bf_edge *edge = find_bf_edge(edges, count, x, y);
    if (!edge) return false;
    if (edge->face[0] != face) {
        *opposite = edge->opposite[0];
        return true;
    }
    if (edge->face[1] != UINT32_MAX) {
        *opposite = edge->opposite[1];
        return true;
    }
    return false;
}

/* Zorin's extraordinary-vertex edge mask. ring is ordered cyclically. */
static geomlib_vec3 extraordinary_edge_point(const geomlib_trimesh *mesh,
                                             const uint32_t *ring, size_t valence,
                                             size_t start, geomlib_vec3 vertex) {
    geomlib_vec3d sum = geomlib_vec3d_mul((geomlib_vec3d){vertex.x, vertex.y, vertex.z}, 0.75);
    for (size_t i = 0; i < valence; i++) {
        size_t index = (start + i) % valence;
        double weight;
        if (valence == 3) {
            static const double table[3] = {5.0 / 12.0, -1.0 / 12.0, -1.0 / 12.0};
            weight = table[i];
        } else if (valence == 4) {
            static const double table[4] = {3.0 / 8.0, 0.0, -1.0 / 8.0, 0.0};
            weight = table[i];
        } else {
            weight = (1.0 / (4.0 * (double)valence)) *
                     (1.0 + 2.0 * cos(2.0 * M_PI * (double)i / (double)valence));
        }
        geomlib_vec3 q = mesh->vertices[ring[index]];
        sum = geomlib_vec3d_add(sum, geomlib_vec3d_mul((geomlib_vec3d){q.x, q.y, q.z}, weight));
    }
    return (geomlib_vec3){(float)sum.x, (float)sum.y, (float)sum.z};
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
    memset(map->buckets, 0, sizeof(map->buckets));
}

static void edge_map_put(edge_map *map, edge_key key, uint32_t new_vertex) {
    uint32_t hash = edge_hash(key);
    edge_map_entry *entry = malloc(sizeof(edge_map_entry));
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
            free(entry);
            entry = next;
        }
    }
}

/* Loop subdivision (one level). Boundary vertices use the cubic B-spline rule
 * 3/4 v + 1/8 (n1 + n2) and boundary edges use their midpoint. When
 * freeze_boundary is true the boundary vertices are held in place instead. */
static int subdivide_loop_internal(
    const geomlib_trimesh *input,
    geomlib_trimesh *output,
    bool freeze_boundary
) {
    if (!input || !output) return GEOMLIB_MESH_ERROR_INVALID;

    /* Build half-edge structure for neighbor queries */
    geomlib_halfedge_mesh hemesh;
    int result = geomlib_halfedge_build(&hemesh, input);
    if (result != GEOMLIB_MESH_SUCCESS) return result;

    size_t edge_count = 0;
    bf_edge *boundary_edges = build_bf_edges(input, &edge_count);
    bool *is_boundary = calloc(input->vertex_count, sizeof(bool));
    uint32_t *boundary_neighbor = malloc(input->vertex_count * 2 * sizeof(uint32_t));
    uint8_t *boundary_count = calloc(input->vertex_count, sizeof(uint8_t));
    if (!boundary_edges || (input->vertex_count && (!is_boundary || !boundary_neighbor ||
                                                   !boundary_count))) {
        free(boundary_edges);
        free(is_boundary);
        free(boundary_neighbor);
        free(boundary_count);
        geomlib_halfedge_free(&hemesh);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (size_t e = 0; e < edge_count; e++) {
        if (boundary_edges[e].face[1] != UINT32_MAX) continue;
        uint32_t a = boundary_edges[e].a, b = boundary_edges[e].b;
        is_boundary[a] = is_boundary[b] = true;
        if (boundary_count[a] < 2) boundary_neighbor[a * 2 + boundary_count[a]++] = b;
        if (boundary_count[b] < 2) boundary_neighbor[b * 2 + boundary_count[b]++] = a;
    }
    free(boundary_edges);

    /* Calculate output mesh size:
     * - Each original vertex remains
     * - Each edge gets a new vertex
     * - Each triangle splits into 4 triangles
     */
    size_t new_vertex_count = input->vertex_count + (hemesh.halfedge_count / 2);
    size_t new_triangle_count = input->triangle_count * 4;
    geomlib_trimesh_init(output, new_vertex_count, new_triangle_count);

    /* Step 1: updated control vertices */
    for (size_t v = 0; v < input->vertex_count; v++) {
        geomlib_vec3 new_pos = input->vertices[v];
        if (is_boundary[v]) {
            if (!freeze_boundary && boundary_count[v] == 2) {
                geomlib_vec3 neighbors = geomlib_vec3_add(
                    input->vertices[boundary_neighbor[v * 2]],
                    input->vertices[boundary_neighbor[v * 2 + 1]]);
                new_pos = geomlib_vec3_add(geomlib_vec3_mul(input->vertices[v], 0.75f),
                                           geomlib_vec3_mul(neighbors, 0.125f));
            }
        } else {
            int valence = 0;
            geomlib_vec3 neighbor_sum = {0.0f, 0.0f, 0.0f};
            for (size_t he = 0; he < hemesh.halfedge_count; he++) {
                if (hemesh.halfedges[he].vertex != v) continue;
                uint32_t next = hemesh.halfedges[he].next;
                if (next == UINT32_MAX) continue;
                neighbor_sum = geomlib_vec3_add(neighbor_sum,
                    input->vertices[hemesh.halfedges[next].vertex]);
                valence++;
            }
            if (valence > 0) {
                float beta = geomlib_loop_beta(valence);
                float alpha = 1.0f - valence * beta;
                new_pos = geomlib_vec3_add(geomlib_vec3_mul(input->vertices[v], alpha),
                                           geomlib_vec3_mul(neighbor_sum, beta));
            }
        }
        geomlib_trimesh_add_vertex(output, new_pos, NULL, NULL);
    }

    /* Step 2: edge vertices */
    edge_map emap;
    edge_map_init(&emap);
    for (size_t t = 0; t < input->triangle_count; t++) {
        uint32_t i0 = input->indices[t * 3 + 0];
        uint32_t i1 = input->indices[t * 3 + 1];
        uint32_t i2 = input->indices[t * 3 + 2];
        edge_key edges[3] = {
            make_edge_key(i0, i1),
            make_edge_key(i1, i2),
            make_edge_key(i2, i0)
        };
        for (int e = 0; e < 3; e++) {
            uint32_t new_v;
            if (edge_map_get(&emap, edges[e], &new_v)) continue;
            geomlib_vec3 v0 = input->vertices[edges[e].v0];
            geomlib_vec3 v1 = input->vertices[edges[e].v1];
            uint32_t opp0 = UINT32_MAX, opp1 = UINT32_MAX;
            for (size_t tt = 0; tt < input->triangle_count; tt++) {
                uint32_t ti[3] = {input->indices[tt * 3], input->indices[tt * 3 + 1],
                                  input->indices[tt * 3 + 2]};
                bool has0 = ti[0] == edges[e].v0 || ti[1] == edges[e].v0 || ti[2] == edges[e].v0;
                bool has1 = ti[0] == edges[e].v1 || ti[1] == edges[e].v1 || ti[2] == edges[e].v1;
                if (!has0 || !has1) continue;
                for (int j = 0; j < 3; j++) {
                    if (ti[j] == edges[e].v0 || ti[j] == edges[e].v1) continue;
                    if (opp0 == UINT32_MAX) opp0 = ti[j];
                    else opp1 = ti[j];
                }
            }
            geomlib_vec3 edge_pos;
            if (opp1 == UINT32_MAX) {
                edge_pos = geomlib_vec3_mul(geomlib_vec3_add(v0, v1), 0.5f);
            } else {
                edge_pos = geomlib_vec3_mul(geomlib_vec3_add(v0, v1), 3.0f / 8.0f);
                edge_pos = geomlib_vec3_add(edge_pos,
                    geomlib_vec3_mul(input->vertices[opp0], 1.0f / 8.0f));
                edge_pos = geomlib_vec3_add(edge_pos,
                    geomlib_vec3_mul(input->vertices[opp1], 1.0f / 8.0f));
            }
            new_v = output->vertex_count;
            geomlib_trimesh_add_vertex(output, edge_pos, NULL, NULL);
            edge_map_put(&emap, edges[e], new_v);
        }
    }

    /* Step 3: split each triangle into four */
    for (size_t t = 0; t < input->triangle_count; t++) {
        uint32_t v0 = input->indices[t * 3 + 0];
        uint32_t v1 = input->indices[t * 3 + 1];
        uint32_t v2 = input->indices[t * 3 + 2];
        /* All triangle edges were inserted into the map in step 2. Initialize
         * defensively so strict builds do not report maybe-uninitialized. */
        uint32_t e01 = 0, e12 = 0, e20 = 0;
        edge_map_get(&emap, make_edge_key(v0, v1), &e01);
        edge_map_get(&emap, make_edge_key(v1, v2), &e12);
        edge_map_get(&emap, make_edge_key(v2, v0), &e20);
        geomlib_trimesh_add_triangle(output, v0, e01, e20);
        geomlib_trimesh_add_triangle(output, v1, e12, e01);
        geomlib_trimesh_add_triangle(output, v2, e20, e12);
        geomlib_trimesh_add_triangle(output, e01, e12, e20);
    }

    edge_map_free(&emap);
    free(is_boundary);
    free(boundary_neighbor);
    free(boundary_count);
    geomlib_halfedge_free(&hemesh);
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_subdivide_loop(
    const geomlib_trimesh *input,
    geomlib_trimesh *output
) {
    return subdivide_loop_internal(input, output, false);
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

/* Modified Butterfly subdivision (Zorin et al.), interpolating. */
int geomlib_subdivide_butterfly(
    const geomlib_trimesh *input,
    geomlib_trimesh *output
) {
    if (!input || !output) return GEOMLIB_MESH_ERROR_INVALID;
    geomlib_vec3 *normals = area_normals(input);
    size_t *offsets = NULL, *counts = NULL;
    uint32_t *ring = NULL;
    size_t edge_count = 0;
    bf_edge *edges = NULL;
    int result = GEOMLIB_MESH_ERROR_MEMORY;
    if (!normals ||
        build_ordered_rings(input, normals, &offsets, &ring, &counts) != GEOMLIB_MESH_SUCCESS)
        goto finish;
    edges = build_bf_edges(input, &edge_count);
    if (!edges) goto finish;
    if (edge_count > UINT32_MAX - input->vertex_count ||
        input->triangle_count > SIZE_MAX / 4) {
        result = GEOMLIB_MESH_ERROR_MEMORY;
        goto finish;
    }
    geomlib_trimesh_init(output, input->vertex_count + edge_count, input->triangle_count * 4);
    if ((input->vertex_count + edge_count && !output->vertices) ||
        (input->triangle_count && !output->indices)) {
        result = GEOMLIB_MESH_ERROR_MEMORY;
        goto finish;
    }

    for (size_t v = 0; v < input->vertex_count; v++) {
        if (geomlib_trimesh_add_vertex(output, input->vertices[v], NULL, NULL) !=
            GEOMLIB_MESH_SUCCESS) goto finish;
    }

    for (size_t e = 0; e < edge_count; e++) {
        bf_edge *edge = &edges[e];
        uint32_t a = edge->a, b = edge->b;
        geomlib_vec3 va = input->vertices[a], vb = input->vertices[b];
        geomlib_vec3 position = geomlib_vec3_mul(geomlib_vec3_add(va, vb), 0.5f);
        if (edge->face[1] != UINT32_MAX) {
            size_t na = counts[a], nb = counts[b];
            if (na == 6 && nb == 6) {
                uint32_t c = edge->opposite[0], d = edge->opposite[1];
                uint32_t wings[4];
                bool full = edge_other_opposite(edges, edge_count, a, c, edge->face[0], &wings[0]) &&
                            edge_other_opposite(edges, edge_count, b, c, edge->face[0], &wings[1]) &&
                            edge_other_opposite(edges, edge_count, a, d, edge->face[1], &wings[2]) &&
                            edge_other_opposite(edges, edge_count, b, d, edge->face[1], &wings[3]);
                position = geomlib_vec3_mul(geomlib_vec3_add(va, vb), 0.5f);
                position = geomlib_vec3_add(position,
                    geomlib_vec3_mul(geomlib_vec3_add(input->vertices[c], input->vertices[d]),
                                     1.0f / 8.0f));
                if (full) {
                    geomlib_vec3d wings_sum = {0};
                    for (int i = 0; i < 4; i++)
                        wings_sum = geomlib_vec3d_add(wings_sum,
                            (geomlib_vec3d){input->vertices[wings[i]].x,
                                            input->vertices[wings[i]].y,
                                            input->vertices[wings[i]].z});
                    position = geomlib_vec3_sub(position, (geomlib_vec3){
                        (float)(wings_sum.x / 16.0), (float)(wings_sum.y / 16.0),
                        (float)(wings_sum.z / 16.0)});
                }
            } else {
                geomlib_vec3d sum = {0};
                int masks = 0;
                if (na != 6 && na > 0) {
                    size_t start = 0;
                    for (size_t i = 0; i < na; i++)
                        if (ring[offsets[a] + i] == b) { start = i; break; }
                    geomlib_vec3 point = extraordinary_edge_point(input, ring + offsets[a], na,
                                                                  start, va);
                    sum = geomlib_vec3d_add(sum, (geomlib_vec3d){point.x, point.y, point.z});
                    masks++;
                }
                if (nb != 6 && nb > 0) {
                    size_t start = 0;
                    for (size_t i = 0; i < nb; i++)
                        if (ring[offsets[b] + i] == a) { start = i; break; }
                    geomlib_vec3 point = extraordinary_edge_point(input, ring + offsets[b], nb,
                                                                  start, vb);
                    sum = geomlib_vec3d_add(sum, (geomlib_vec3d){point.x, point.y, point.z});
                    masks++;
                }
                if (masks) position = (geomlib_vec3){(float)(sum.x / masks),
                                                     (float)(sum.y / masks),
                                                     (float)(sum.z / masks)};
            }
        }
        if (geomlib_trimesh_add_vertex(output, position, NULL, NULL) != GEOMLIB_MESH_SUCCESS)
            goto finish;
        edge->new_vertex = (uint32_t)(output->vertex_count - 1);
    }

    for (size_t t = 0; t < input->triangle_count; t++) {
        uint32_t v0 = input->indices[t * 3], v1 = input->indices[t * 3 + 1],
                 v2 = input->indices[t * 3 + 2];
        bf_edge *e01 = find_bf_edge(edges, edge_count, v0, v1);
        bf_edge *e12 = find_bf_edge(edges, edge_count, v1, v2);
        bf_edge *e20 = find_bf_edge(edges, edge_count, v2, v0);
        if (!e01 || !e12 || !e20) continue;
        uint32_t m01 = e01->new_vertex, m12 = e12->new_vertex, m20 = e20->new_vertex;
        int status = geomlib_trimesh_add_triangle(output, v0, m01, m20);
        if (status == GEOMLIB_MESH_SUCCESS) status = geomlib_trimesh_add_triangle(output, v1, m12, m01);
        if (status == GEOMLIB_MESH_SUCCESS) status = geomlib_trimesh_add_triangle(output, v2, m20, m12);
        if (status == GEOMLIB_MESH_SUCCESS) status = geomlib_trimesh_add_triangle(output, m01, m12, m20);
        if (status != GEOMLIB_MESH_SUCCESS) { result = status; goto finish; }
    }

    result = GEOMLIB_MESH_SUCCESS;

finish:
    if (result != GEOMLIB_MESH_SUCCESS) geomlib_trimesh_free(output);
    free(normals);
    free(offsets);
    free(ring);
    free(counts);
    free(edges);
    return result;
}

/* Project every interior vertex onto the Loop limit surface in one pass. */
static bool boundary_vertex_from_edges(const bf_edge *edges, size_t count, uint32_t vertex);
static void loop_project_to_limit(geomlib_trimesh *mesh) {
    if (!mesh || !mesh->vertex_count) return;
    geomlib_vec3 *normals = area_normals(mesh);
    size_t edge_count = 0;
    bf_edge *edges = build_bf_edges(mesh, &edge_count);
    size_t *offsets = NULL, *counts = NULL;
    uint32_t *ring = NULL;
    geomlib_vec3 *positions = malloc(mesh->vertex_count * sizeof(*positions));
    bool *is_boundary = calloc(mesh->vertex_count, sizeof(*is_boundary));
    if (!normals || !edges || !positions || !is_boundary ||
        build_ordered_rings(mesh, normals, &offsets, &ring, &counts) != GEOMLIB_MESH_SUCCESS) {
        free(positions);
        free(is_boundary);
        free(normals);
        free(edges);
        free(offsets);
        free(ring);
        free(counts);
        return;
    }
    for (size_t e = 0; e < edge_count; e++) {
        if (edges[e].face[1] != UINT32_MAX) continue;
        is_boundary[edges[e].a] = true;
        is_boundary[edges[e].b] = true;
    }
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        positions[v] = mesh->vertices[v];
        size_t valence = counts[v];
        if (valence < 3 || is_boundary[v]) continue;
        double alpha = pow(3.0 / 8.0 + cos(2.0 * M_PI / (double)valence) / 4.0, 2.0) + 3.0 / 8.0;
        double epsilon = 3.0 / (11.0 - 8.0 * alpha);
        geomlib_vec3d mean = {0};
        for (size_t i = 0; i < valence; i++) {
            geomlib_vec3 p = mesh->vertices[ring[offsets[v] + i]];
            mean = geomlib_vec3d_add(mean, (geomlib_vec3d){p.x, p.y, p.z});
        }
        mean = geomlib_vec3d_mul(mean, 1.0 / (double)valence);
        geomlib_vec3 current = mesh->vertices[v];
        geomlib_vec3d limit = geomlib_vec3d_add(
            geomlib_vec3d_mul((geomlib_vec3d){current.x, current.y, current.z}, (float)epsilon),
            geomlib_vec3d_mul(mean, (float)(1.0 - epsilon)));
        positions[v] = (geomlib_vec3){(float)limit.x, (float)limit.y, (float)limit.z};
    }
    memcpy(mesh->vertices, positions, mesh->vertex_count * sizeof(*positions));
    free(positions);
    free(is_boundary);
    free(normals);
    free(edges);
    free(offsets);
    free(ring);
    free(counts);
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
                result = subdivide_loop_internal(current, next, params->preserve_boundaries);
                break;
            case GEOMLIB_SUBDIV_BUTTERFLY:
                result = geomlib_subdivide_butterfly(current, next);
                break;
            case GEOMLIB_SUBDIV_CATMULL_CLARK:
                result = geomlib_subdivide_catmull_clark(current, next);
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

    if (params->scheme == GEOMLIB_SUBDIV_LOOP && params->project_to_limit)
        loop_project_to_limit(current);

    /* Copy result to output */
    *output = *current;

    return GEOMLIB_MESH_SUCCESS;
}

static bool boundary_vertex_from_edges(const bf_edge *edges, size_t count, uint32_t vertex) {
    for (size_t i = 0; i < count; i++)
        if (edges[i].face[1] == UINT32_MAX &&
            (edges[i].a == vertex || edges[i].b == vertex)) return true;
    return false;
}

/* Loop limit position via the evaluation mask
 * alpha_n = (3/8 + cos(2*pi/n)/4)^2 + 3/8, epsilon_n = 3/(11 - 8*alpha_n). */
geomlib_vec3 geomlib_loop_limit_position(
    const geomlib_trimesh *mesh,
    uint32_t vertex_idx
) {
    if (!mesh || vertex_idx >= mesh->vertex_count)
        return (geomlib_vec3){0.0f, 0.0f, 0.0f};
    geomlib_vec3 result = mesh->vertices[vertex_idx];
    geomlib_vec3 *normals = area_normals(mesh);
    size_t edge_count = 0;
    bf_edge *edges = build_bf_edges(mesh, &edge_count);
    size_t *offsets = NULL, *counts = NULL;
    uint32_t *ring = NULL;
    if (!normals || !edges ||
        build_ordered_rings(mesh, normals, &offsets, &ring, &counts) != GEOMLIB_MESH_SUCCESS)
        goto finish;
    size_t valence = counts[vertex_idx];
    if (valence < 3 || boundary_vertex_from_edges(edges, edge_count, vertex_idx)) goto finish;

    double alpha = pow(3.0 / 8.0 + cos(2.0 * M_PI / (double)valence) / 4.0, 2.0) + 3.0 / 8.0;
    double epsilon = 3.0 / (11.0 - 8.0 * alpha);
    geomlib_vec3d mean = {0};
    for (size_t i = 0; i < valence; i++) {
        geomlib_vec3 p = mesh->vertices[ring[offsets[vertex_idx] + i]];
        mean = geomlib_vec3d_add(mean, (geomlib_vec3d){p.x, p.y, p.z});
    }
    mean = geomlib_vec3d_mul(mean, 1.0 / (double)valence);
    geomlib_vec3d v = {result.x, result.y, result.z};
    geomlib_vec3d limit = geomlib_vec3d_add(geomlib_vec3d_mul(v, (float)epsilon),
                                           geomlib_vec3d_mul(mean, (float)(1.0 - epsilon)));
    result = (geomlib_vec3){(float)limit.x, (float)limit.y, (float)limit.z};

finish:
    free(normals);
    free(edges);
    free(offsets);
    free(ring);
    free(counts);
    return result;
}

/* Loop limit normal from the tangent masks tau_k(1)=cos(2*pi*k/n),
 * tau_k(2)=sin(2*pi*k/n) over the cyclically ordered one-ring. */
geomlib_vec3 geomlib_loop_limit_normal(
    const geomlib_trimesh *mesh,
    uint32_t vertex_idx
) {
    if (!mesh || vertex_idx >= mesh->vertex_count)
        return (geomlib_vec3){0.0f, 1.0f, 0.0f};
    geomlib_vec3 fallback = mesh->normals && vertex_idx < mesh->vertex_count
        ? mesh->normals[vertex_idx] : geomlib_vec3_create(0, 1, 0);
    geomlib_vec3 *normals = area_normals(mesh);
    size_t edge_count = 0;
    bf_edge *edges = build_bf_edges(mesh, &edge_count);
    size_t *offsets = NULL, *counts = NULL;
    uint32_t *ring = NULL;
    geomlib_vec3 result = fallback;
    if (!normals || !edges ||
        build_ordered_rings(mesh, normals, &offsets, &ring, &counts) != GEOMLIB_MESH_SUCCESS)
        goto finish;
    size_t valence = counts[vertex_idx];
    if (valence < 3 || boundary_vertex_from_edges(edges, edge_count, vertex_idx)) {
        result = normals[vertex_idx];
        goto finish;
    }
    geomlib_vec3d tangent1 = {0}, tangent2 = {0};
    for (size_t k = 0; k < valence; k++) {
        geomlib_vec3 p = mesh->vertices[ring[offsets[vertex_idx] + k]];
        double angle = 2.0 * M_PI * (double)k / (double)valence;
        tangent1 = geomlib_vec3d_add(tangent1,
            geomlib_vec3d_mul((geomlib_vec3d){p.x, p.y, p.z}, cos(angle)));
        tangent2 = geomlib_vec3d_add(tangent2,
            geomlib_vec3d_mul((geomlib_vec3d){p.x, p.y, p.z}, sin(angle)));
    }
    geomlib_vec3 n = geomlib_vec3_cross((geomlib_vec3){(float)tangent1.x, (float)tangent1.y,
                                                       (float)tangent1.z},
                                        (geomlib_vec3){(float)tangent2.x, (float)tangent2.y,
                                                       (float)tangent2.z});
    float length = geomlib_vec3_length(n);
    if (!(length > 1e-12f)) {
        result = normals[vertex_idx];
        goto finish;
    }
    n = geomlib_vec3_mul(n, 1.0f / length);
    result = geomlib_vec3_dot(n, normals[vertex_idx]) < 0.0f ? geomlib_vec3_mul(n, -1.0f) : n;

finish:
    free(normals);
    free(edges);
    free(offsets);
    free(ring);
    free(counts);
    return result;
}
