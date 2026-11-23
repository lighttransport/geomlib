/*
 * geomlib - Mesh Topology Analysis Implementation
 */

#include "geomlib/geomlib_topology.h"
#include "geomlib/geomlib_bvh.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

/* ============================================================================
 * Edge Hash Map for Fast Edge Lookup
 * ============================================================================ */

typedef struct edge_key {
    uint32_t v0, v1;  /* v0 < v1 */
} edge_key;

typedef struct edge_hash_entry {
    edge_key key;
    uint32_t edge_index;
    uint32_t face_count;
    uint32_t forward_count;
    uint32_t faces[2];
    uint32_t first_corner[2];
    struct edge_hash_entry *next;
} edge_hash_entry;

typedef struct edge_hash_map {
    edge_hash_entry **buckets;
    size_t bucket_count;
    size_t entry_count;
} edge_hash_map;

static uint32_t hash_edge(edge_key key, size_t bucket_count) {
    uint32_t h = 2166136261u;
    h = (h ^ key.v0) * 16777619u;
    h = (h ^ key.v1) * 16777619u;
    return h % bucket_count;
}

static edge_hash_map edge_hash_map_create(size_t capacity) {
    edge_hash_map map;
    map.bucket_count = capacity > 0 ? capacity : 1024;
    map.buckets = (edge_hash_entry**)calloc(map.bucket_count, sizeof(edge_hash_entry*));
    map.entry_count = 0;
    return map;
}

static void edge_hash_map_free(edge_hash_map *map) {
    for (size_t i = 0; map->buckets && i < map->bucket_count; i++) {
        edge_hash_entry *entry = map->buckets[i];
        while (entry) {
            edge_hash_entry *next = entry->next;
            free(entry);
            entry = next;
        }
    }
    free(map->buckets);
}

static edge_hash_entry* edge_hash_map_find(edge_hash_map *map, edge_key key) {
    uint32_t bucket = hash_edge(key, map->bucket_count);
    edge_hash_entry *entry = map->buckets[bucket];

    while (entry) {
        if (entry->key.v0 == key.v0 && entry->key.v1 == key.v1) {
            return entry;
        }
        entry = entry->next;
    }
    return NULL;
}

static edge_hash_entry* edge_hash_map_insert(edge_hash_map *map, edge_key key) {
    edge_hash_entry *existing = edge_hash_map_find(map, key);
    if (existing) {
        return existing;
    }

    uint32_t bucket = hash_edge(key, map->bucket_count);
    edge_hash_entry *entry = (edge_hash_entry*)calloc(1, sizeof(edge_hash_entry));
    if (!entry) return NULL;
    entry->key = key;
    entry->edge_index = (uint32_t)map->entry_count++;
    entry->face_count = 0;
    entry->next = map->buckets[bucket];
    map->buckets[bucket] = entry;

    return entry;
}

static edge_key make_edge_key(uint32_t v0, uint32_t v1) {
    edge_key key;
    if (v0 < v1) {
        key.v0 = v0;
        key.v1 = v1;
    } else {
        key.v0 = v1;
        key.v1 = v0;
    }
    return key;
}

/* ============================================================================
 * Edge Mesh Construction
 * ============================================================================ */

static uint32_t find_root(uint32_t *parent, uint32_t x);
static void union_sets(uint32_t *parent, uint32_t *rank, uint32_t x, uint32_t y);

static bool build_edge_map(const geomlib_trimesh *mesh, edge_hash_map *map,
                           uint32_t *parent, uint32_t *rank) {
    if (!mesh || mesh->triangle_count > UINT32_MAX / 3 ||
        (mesh->triangle_count && !mesh->indices)) return false;
    *map = edge_hash_map_create(mesh->triangle_count * 3);
    if (!map->buckets) return false;

    for (size_t face = 0; face < mesh->triangle_count; face++) {
        for (int i = 0; i < 3; i++) {
            uint32_t a = (uint32_t)(face * 3 + i);
            uint32_t b = (uint32_t)(face * 3 + (i + 1) % 3);
            uint32_t va = mesh->indices[a], vb = mesh->indices[b];
            if (va >= mesh->vertex_count || vb >= mesh->vertex_count) return false;
            edge_hash_entry *edge = edge_hash_map_insert(map, make_edge_key(va, vb));
            if (!edge) return false;

            bool forward = va < vb;
            uint32_t low = forward ? a : b;
            uint32_t high = forward ? b : a;
            if (!edge->face_count) {
                edge->first_corner[0] = low;
                edge->first_corner[1] = high;
            } else if (parent) {
                union_sets(parent, rank, low, edge->first_corner[0]);
                union_sets(parent, rank, high, edge->first_corner[1]);
            }
            if (edge->face_count < 2) edge->faces[edge->face_count] = (uint32_t)face;
            edge->face_count++;
            edge->forward_count += forward;
        }
    }
    return true;
}

geomlib_edge_mesh geomlib_build_edge_mesh(const geomlib_trimesh *mesh) {
    geomlib_edge_mesh result = {0};
    edge_hash_map map = {0};
    if (!build_edge_map(mesh, &map, NULL, NULL)) {
        edge_hash_map_free(&map);
        return result;
    }

    result.edge_count = map.entry_count;
    result.vertex_count = mesh->vertex_count;
    result.edges = (geomlib_edge*)malloc(result.edge_count * sizeof(geomlib_edge));
    result.vertex_edges = (uint32_t**)calloc(result.vertex_count, sizeof(uint32_t*));
    result.vertex_edge_count = (uint32_t*)calloc(result.vertex_count, sizeof(uint32_t));
    if ((result.edge_count && !result.edges) ||
        (result.vertex_count && (!result.vertex_edges || !result.vertex_edge_count))) goto error;

    size_t index = 0;
    for (size_t i = 0; i < map.bucket_count; i++) {
        for (edge_hash_entry *edge = map.buckets[i]; edge; edge = edge->next) {
            geomlib_edge *out = &result.edges[index++];
            out->v0 = edge->key.v0;
            out->v1 = edge->key.v1;
            out->face0 = edge->faces[0];
            out->face1 = edge->face_count > 1 ? edge->faces[1] : UINT32_MAX;
            result.vertex_edge_count[out->v0]++;
            result.vertex_edge_count[out->v1]++;
        }
    }

    for (size_t i = 0; i < result.vertex_count; i++) {
        if (!result.vertex_edge_count[i]) continue;
        result.vertex_edges[i] = (uint32_t*)malloc(result.vertex_edge_count[i] * sizeof(uint32_t));
        if (!result.vertex_edges[i]) goto error;
        result.vertex_edge_count[i] = 0;
    }
    for (size_t i = 0; i < result.edge_count; i++) {
        uint32_t a = result.edges[i].v0, b = result.edges[i].v1;
        result.vertex_edges[a][result.vertex_edge_count[a]++] = (uint32_t)i;
        result.vertex_edges[b][result.vertex_edge_count[b]++] = (uint32_t)i;
    }
    edge_hash_map_free(&map);
    return result;

error:
    edge_hash_map_free(&map);
    geomlib_edge_mesh_free(&result);
    return result;
}

void geomlib_edge_mesh_free(geomlib_edge_mesh *edge_mesh) {
    free(edge_mesh->edges);
    for (size_t i = 0; edge_mesh->vertex_edges && i < edge_mesh->vertex_count; i++) {
        free(edge_mesh->vertex_edges[i]);
    }
    free(edge_mesh->vertex_edges);
    free(edge_mesh->vertex_edge_count);
    memset(edge_mesh, 0, sizeof(*edge_mesh));
}

/* ============================================================================
 * Euler Characteristic
 * ============================================================================ */

geomlib_euler_characteristic geomlib_compute_euler_characteristic(
    const geomlib_trimesh *mesh
) {
    geomlib_euler_characteristic result = {0};
    if (!mesh) return result;

    result.vertices = (int32_t)mesh->vertex_count;
    result.faces = (int32_t)mesh->triangle_count;

    /* Build edge mesh to count edges */
    geomlib_edge_mesh edge_mesh = geomlib_build_edge_mesh(mesh);
    result.edges = (int32_t)edge_mesh.edge_count;

    /* Compute connected components */
    geomlib_component_info components = geomlib_find_connected_components(mesh);
    result.components = components.component_count;
    geomlib_component_info_free(&components);

    /* Extract boundary loops */
    geomlib_boundary_loops loops = geomlib_extract_boundary_loops(mesh);
    result.boundaries = (int32_t)loops.loop_count;
    geomlib_boundary_loops_free(&loops);

    /* Euler characteristic: χ = V - E + F */
    result.chi = result.vertices - result.edges + result.faces;

    /* Genus for orientable surface: g = (2 - χ - b) / 2 */
    /* Adjusted for multiple components: χ = 2c - 2g - b */
    result.genus = (2 * result.components - result.chi - result.boundaries) / 2;

    geomlib_edge_mesh_free(&edge_mesh);

    return result;
}

/* ============================================================================
 * Boundary Detection
 * ============================================================================ */

bool geomlib_is_boundary_vertex(
    const geomlib_trimesh *mesh,
    uint32_t vertex_idx
) {
    geomlib_edge_mesh edge_mesh = geomlib_build_edge_mesh(mesh);

    bool is_boundary = false;
    if (vertex_idx >= edge_mesh.vertex_count) {
        geomlib_edge_mesh_free(&edge_mesh);
        return false;
    }
    for (uint32_t i = 0; i < edge_mesh.vertex_edge_count[vertex_idx]; i++) {
        uint32_t edge_idx = edge_mesh.vertex_edges[vertex_idx][i];
        if (edge_mesh.edges[edge_idx].face1 == UINT32_MAX) {
            is_boundary = true;
            break;
        }
    }

    geomlib_edge_mesh_free(&edge_mesh);
    return is_boundary;
}

geomlib_boundary_loops geomlib_extract_boundary_loops(const geomlib_trimesh *mesh) {
    geomlib_boundary_loops result = {0};
    geomlib_edge_mesh edges = geomlib_build_edge_mesh(mesh);
    if (!edges.edge_count || !mesh->vertices) {
        geomlib_edge_mesh_free(&edges);
        return result;
    }
    uint32_t *degree = (uint32_t*)calloc(edges.vertex_count, sizeof(uint32_t));
    bool *visited = (bool*)calloc(edges.edge_count, sizeof(bool));
    uint32_t *path = (uint32_t*)malloc(edges.edge_count * sizeof(uint32_t));
    if (!degree || !visited || !path) goto finish;

    for (size_t i = 0; i < edges.edge_count; i++) {
        if (edges.edges[i].face1 != UINT32_MAX) continue;
        degree[edges.edges[i].v0]++;
        degree[edges.edges[i].v1]++;
    }

    size_t capacity = 0;
    for (size_t i = 0; i < edges.edge_count; i++) {
        if (visited[i] || edges.edges[i].face1 != UINT32_MAX) continue;
        uint32_t start = edges.edges[i].v0, vertex = start;
        if (degree[start] != 2) continue;
        uint32_t edge = (uint32_t)i;
        size_t count = 0;
        float length = 0.0f;
        do {
            visited[edge] = true;
            path[count++] = vertex;
            geomlib_edge current = edges.edges[edge];
            uint32_t next = current.v0 == vertex ? current.v1 : current.v0;
            length += geomlib_vec3_length(geomlib_vec3_sub(mesh->vertices[vertex], mesh->vertices[next]));
            vertex = next;
            if (vertex == start || degree[vertex] != 2) break;
            edge = UINT32_MAX;
            for (uint32_t j = 0; j < edges.vertex_edge_count[vertex]; j++) {
                uint32_t candidate = edges.vertex_edges[vertex][j];
                if (!visited[candidate] && edges.edges[candidate].face1 == UINT32_MAX) {
                    edge = candidate;
                    break;
                }
            }
        } while (edge != UINT32_MAX);

        if (vertex != start || count < 3) continue;
        if (result.loop_count == capacity) {
            size_t next_capacity = capacity ? capacity * 2 : 4;
            geomlib_boundary_loop *loops = (geomlib_boundary_loop*)realloc(
                result.loops, next_capacity * sizeof(*loops));
            if (!loops) goto error;
            result.loops = loops;
            capacity = next_capacity;
        }
        uint32_t *vertices = (uint32_t*)malloc(count * sizeof(uint32_t));
        if (!vertices) goto error;
        memcpy(vertices, path, count * sizeof(uint32_t));
        result.loops[result.loop_count++] = (geomlib_boundary_loop){vertices, count, length};
    }
    goto finish;
error:
    geomlib_boundary_loops_free(&result);
    result = (geomlib_boundary_loops){0};
finish:
    free(degree);
    free(visited);
    free(path);
    geomlib_edge_mesh_free(&edges);
    return result;
}

void geomlib_boundary_loops_free(geomlib_boundary_loops *loops) {
    for (size_t i = 0; i < loops->loop_count; i++) {
        free(loops->loops[i].vertices);
    }
    free(loops->loops);
}

/* ============================================================================
 * Connected Components (Union-Find)
 * ============================================================================ */

static uint32_t find_root(uint32_t *parent, uint32_t x) {
    if (parent[x] != x) {
        parent[x] = find_root(parent, parent[x]);  /* Path compression */
    }
    return parent[x];
}

static void union_sets(uint32_t *parent, uint32_t *rank, uint32_t x, uint32_t y) {
    uint32_t root_x = find_root(parent, x);
    uint32_t root_y = find_root(parent, y);

    if (root_x != root_y) {
        /* Union by rank */
        if (rank[root_x] < rank[root_y]) {
            parent[root_x] = root_y;
        } else if (rank[root_x] > rank[root_y]) {
            parent[root_y] = root_x;
        } else {
            parent[root_y] = root_x;
            rank[root_x]++;
        }
    }
}

geomlib_component_info geomlib_find_connected_components(
    const geomlib_trimesh *mesh
) {
    geomlib_component_info info;

    /* Initialize union-find */
    uint32_t *parent = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    uint32_t *rank = (uint32_t*)calloc(mesh->vertex_count, sizeof(uint32_t));

    for (size_t i = 0; i < mesh->vertex_count; i++) {
        parent[i] = (uint32_t)i;
    }

    /* Union vertices in each face */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t v0 = mesh->indices[i * 3 + 0];
        uint32_t v1 = mesh->indices[i * 3 + 1];
        uint32_t v2 = mesh->indices[i * 3 + 2];

        union_sets(parent, rank, v0, v1);
        union_sets(parent, rank, v1, v2);
    }

    /* Find roots and assign component IDs */
    uint32_t *root_to_component = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        root_to_component[i] = UINT32_MAX;
    }

    info.vertex_labels = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    info.component_count = 0;

    for (size_t i = 0; i < mesh->vertex_count; i++) {
        uint32_t root = find_root(parent, (uint32_t)i);
        if (root_to_component[root] == UINT32_MAX) {
            root_to_component[root] = info.component_count++;
        }
        info.vertex_labels[i] = root_to_component[root];
    }

    /* Assign face labels */
    info.face_labels = (uint32_t*)malloc(mesh->triangle_count * sizeof(uint32_t));
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        info.face_labels[i] = info.vertex_labels[mesh->indices[i * 3]];
    }

    /* Count vertices and faces per component */
    info.component_vertex_counts = (uint32_t*)calloc(info.component_count, sizeof(uint32_t));
    info.component_face_counts = (uint32_t*)calloc(info.component_count, sizeof(uint32_t));

    for (size_t i = 0; i < mesh->vertex_count; i++) {
        info.component_vertex_counts[info.vertex_labels[i]]++;
    }
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        info.component_face_counts[info.face_labels[i]]++;
    }

    free(parent);
    free(rank);
    free(root_to_component);

    return info;
}

void geomlib_component_info_free(geomlib_component_info *info) {
    free(info->vertex_labels);
    free(info->face_labels);
    free(info->component_vertex_counts);
    free(info->component_face_counts);
}

int geomlib_extract_largest_component(
    const geomlib_trimesh *mesh,
    geomlib_trimesh *output
) {
    geomlib_component_info info = geomlib_find_connected_components(mesh);

    /* Find largest component */
    uint32_t largest = 0;
    for (uint32_t i = 1; i < info.component_count; i++) {
        if (info.component_face_counts[i] > info.component_face_counts[largest]) {
            largest = i;
        }
    }

    /* Build vertex remap */
    uint32_t *vertex_remap = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    uint32_t new_vertex_count = 0;

    for (size_t i = 0; i < mesh->vertex_count; i++) {
        if (info.vertex_labels[i] == largest) {
            vertex_remap[i] = new_vertex_count++;
        } else {
            vertex_remap[i] = UINT32_MAX;
        }
    }

    /* Copy vertices */
    output->vertices = (geomlib_vec3*)malloc(new_vertex_count * sizeof(geomlib_vec3));
    output->vertex_count = new_vertex_count;
    output->vertex_capacity = new_vertex_count;

    for (size_t i = 0; i < mesh->vertex_count; i++) {
        if (vertex_remap[i] != UINT32_MAX) {
            output->vertices[vertex_remap[i]] = mesh->vertices[i];
        }
    }

    /* Copy faces */
    uint32_t face_count = info.component_face_counts[largest];
    output->indices = (uint32_t*)malloc(face_count * 3 * sizeof(uint32_t));
    output->triangle_count = face_count;
    output->triangle_capacity = face_count;

    uint32_t face_idx = 0;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        if (info.face_labels[i] == largest) {
            output->indices[face_idx * 3 + 0] = vertex_remap[mesh->indices[i * 3 + 0]];
            output->indices[face_idx * 3 + 1] = vertex_remap[mesh->indices[i * 3 + 1]];
            output->indices[face_idx * 3 + 2] = vertex_remap[mesh->indices[i * 3 + 2]];
            face_idx++;
        }
    }

    output->normals = NULL;
    output->texcoords = NULL;

    free(vertex_remap);
    geomlib_component_info_free(&info);

    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * Manifold Detection
 * ============================================================================ */

/* Corners joined by an edge at a vertex belong to the same incident fan. */
static uint32_t *build_corner_fans(const geomlib_trimesh *mesh, edge_hash_map *map) {
    if (mesh->triangle_count > UINT32_MAX / 3 ||
        mesh->triangle_count > SIZE_MAX / (3 * sizeof(uint32_t))) return NULL;
    size_t count = mesh->triangle_count * 3;
    uint32_t *parent = (uint32_t*)malloc(count * sizeof(uint32_t));
    uint32_t *rank = (uint32_t*)calloc(count, sizeof(uint32_t));
    if (!parent || !rank) {
        free(parent);
        free(rank);
        return NULL;
    }
    for (size_t i = 0; i < count; i++) parent[i] = (uint32_t)i;
    bool ok = build_edge_map(mesh, map, parent, rank);
    free(rank);
    if (!ok) {
        free(parent);
        return NULL;
    }
    for (size_t i = 0; i < count; i++) parent[i] = find_root(parent, (uint32_t)i);
    return parent;
}

geomlib_manifold_analysis geomlib_analyze_manifold(const geomlib_trimesh *mesh) {
    geomlib_manifold_analysis result = {0};
    if (!mesh) return result;
    result.is_manifold = result.is_edge_manifold = result.is_vertex_manifold = true;
    result.is_closed = result.is_oriented = true;
    if (!mesh->triangle_count) return result;

    edge_hash_map map = {0};
    uint32_t *fans = build_corner_fans(mesh, &map);
    uint32_t *first = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    uint32_t *incident = (uint32_t*)calloc(mesh->vertex_count, sizeof(uint32_t));
    bool *nonmanifold = (bool*)calloc(mesh->vertex_count, sizeof(bool));
    if (!fans || !first || !incident || !nonmanifold) goto error;
    memset(first, 0xff, mesh->vertex_count * sizeof(uint32_t));
    for (size_t i = 0; i < mesh->triangle_count * 3; i++) {
        uint32_t vertex = mesh->indices[i];
        if (first[vertex] == UINT32_MAX) first[vertex] = fans[i];
        if (first[vertex] != fans[i]) nonmanifold[vertex] = true;
        incident[vertex]++;
    }

    result.nonmanifold_elements = (geomlib_nonmanifold_element*)malloc(
        (map.entry_count + mesh->vertex_count) * sizeof(geomlib_nonmanifold_element));
    if (!result.nonmanifold_elements) goto error;
    uint32_t index = 0;
    for (size_t i = 0; i < map.bucket_count; i++) {
        for (edge_hash_entry *edge = map.buckets[i]; edge; edge = edge->next, index++) {
            if (edge->face_count == 1) result.is_closed = false;
            if (edge->face_count == 2 && edge->forward_count != 1) result.is_oriented = false;
            if (edge->face_count > 2) {
                result.is_edge_manifold = false;
                result.is_oriented = false;
                nonmanifold[edge->key.v0] = nonmanifold[edge->key.v1] = true;
                result.nonmanifold_elements[result.nonmanifold_count++] =
                    (geomlib_nonmanifold_element){GEOMLIB_NONMANIFOLD_EDGE, index, edge->face_count};
            }
        }
    }
    for (size_t v = 0; v < mesh->vertex_count; v++) {
        if (!nonmanifold[v]) continue;
        result.is_vertex_manifold = false;
        result.nonmanifold_elements[result.nonmanifold_count++] =
            (geomlib_nonmanifold_element){GEOMLIB_NONMANIFOLD_VERTEX, (uint32_t)v, incident[v]};
    }
    result.is_manifold = result.is_edge_manifold && result.is_vertex_manifold;
    if (!result.nonmanifold_count) {
        free(result.nonmanifold_elements);
        result.nonmanifold_elements = NULL;
    }
    goto finish;
error:
    result = (geomlib_manifold_analysis){0};
finish:
    free(fans);
    free(first);
    free(incident);
    free(nonmanifold);
    edge_hash_map_free(&map);
    return result;
}

void geomlib_manifold_analysis_free(geomlib_manifold_analysis *analysis) {
    free(analysis->nonmanifold_elements);
}

/* ============================================================================
 * Vertex Adjacency
 * ============================================================================ */

geomlib_vertex_adjacency geomlib_build_vertex_adjacency(const geomlib_trimesh *mesh) {
    geomlib_vertex_adjacency result = {0};
    if (!mesh || !mesh->vertex_count || mesh->triangle_count > UINT32_MAX / 6 ||
        mesh->triangle_count > SIZE_MAX / (6 * sizeof(uint32_t)) ||
        (mesh->triangle_count && !mesh->indices)) return result;

    size_t count = mesh->triangle_count * 3;
    result.neighbors = (uint32_t**)calloc(mesh->vertex_count, sizeof(uint32_t*));
    result.faces = (uint32_t**)calloc(mesh->vertex_count, sizeof(uint32_t*));
    result.neighbor_counts = (uint32_t*)calloc(mesh->vertex_count, sizeof(uint32_t));
    result.face_counts = (uint32_t*)calloc(mesh->vertex_count, sizeof(uint32_t));
    if (!result.neighbors || !result.faces || !result.neighbor_counts || !result.face_counts) goto error;
    if (!count) return result;

    result.neighbors[0] = (uint32_t*)malloc(count * 2 * sizeof(uint32_t));
    result.faces[0] = (uint32_t*)malloc(count * sizeof(uint32_t));
    if (!result.neighbors[0] || !result.faces[0]) goto error;

    for (size_t i = 0; i < count; i++) {
        uint32_t vertex = mesh->indices[i];
        if (vertex >= mesh->vertex_count) goto error;
        result.face_counts[vertex]++;
    }
    size_t offset = 0;
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        result.neighbors[i] = result.neighbors[0] + offset * 2;
        result.faces[i] = result.faces[0] + offset;
        offset += result.face_counts[i];
        result.face_counts[i] = 0;
    }
    for (size_t face = 0; face < mesh->triangle_count; face++) {
        const uint32_t *vertices = &mesh->indices[face * 3];
        for (int i = 0; i < 3; i++) {
            uint32_t vertex = vertices[i];
            result.neighbors[vertex][result.neighbor_counts[vertex]++] = vertices[(i + 1) % 3];
            result.neighbors[vertex][result.neighbor_counts[vertex]++] = vertices[(i + 2) % 3];
            result.faces[vertex][result.face_counts[vertex]++] = (uint32_t)face;
        }
    }
    return result;

error:
    geomlib_vertex_adjacency_free(&result);
    return result;
}

void geomlib_vertex_adjacency_free(geomlib_vertex_adjacency *adjacency) {
    if (adjacency->neighbors) free(adjacency->neighbors[0]);
    if (adjacency->faces) free(adjacency->faces[0]);
    free(adjacency->neighbors);
    free(adjacency->neighbor_counts);
    free(adjacency->faces);
    free(adjacency->face_counts);
    memset(adjacency, 0, sizeof(*adjacency));
}

/* ============================================================================
 * Valence Statistics
 * ============================================================================ */

geomlib_valence_stats geomlib_compute_valence_stats(
    const geomlib_trimesh *mesh
) {
    geomlib_valence_stats stats = {0};
    if (!mesh || !mesh->vertex_count) return stats;

    geomlib_edge_mesh edges = geomlib_build_edge_mesh(mesh);
    if (!edges.vertex_edge_count) return stats;

    stats.valences = edges.vertex_edge_count;
    edges.vertex_edge_count = NULL;
    geomlib_edge_mesh_free(&edges);

    stats.min_valence = UINT32_MAX;
    size_t total_valence = 0;
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        uint32_t count = stats.valences[i];
        if (count < stats.min_valence) stats.min_valence = count;
        if (count > stats.max_valence) stats.max_valence = count;
        total_valence += count;
    }
    stats.avg_valence = (float)total_valence / (float)mesh->vertex_count;

    stats.histogram = (uint32_t*)calloc((size_t)stats.max_valence + 1, sizeof(uint32_t));
    if (stats.histogram) {
        stats.histogram_size = (size_t)stats.max_valence + 1;
        for (size_t i = 0; i < mesh->vertex_count; i++) stats.histogram[stats.valences[i]]++;
    }

    return stats;
}

void geomlib_valence_stats_free(geomlib_valence_stats *stats) {
    free(stats->valences);
    free(stats->histogram);
}

/* ============================================================================
 * Self-Intersection Detection (Basic BVH-based)
 * ============================================================================ */

/* Helper: Test if two triangles intersect (basic AABB and edge-triangle tests) */
static geomlib_aabb triangle_bounds(const geomlib_vec3 triangle[3]) {
    geomlib_aabb bounds = {triangle[0], triangle[0]};
    for (int i = 1; i < 3; i++) {
        geomlib_vec3 p = triangle[i];
        bounds.min.x = fminf(bounds.min.x, p.x);
        bounds.min.y = fminf(bounds.min.y, p.y);
        bounds.min.z = fminf(bounds.min.z, p.z);
        bounds.max.x = fmaxf(bounds.max.x, p.x);
        bounds.max.y = fmaxf(bounds.max.y, p.y);
        bounds.max.z = fmaxf(bounds.max.z, p.z);
    }
    return bounds;
}

static float triangle_half_plane(geomlib_vec3 a, geomlib_vec3 b,
                                 geomlib_vec3 normal, geomlib_vec3 point) {
    geomlib_vec3 edge = geomlib_vec3_sub(b, a);
    geomlib_vec3 inward = geomlib_vec3_normalize(geomlib_vec3_cross(normal, edge));
    return geomlib_vec3_dot(inward, geomlib_vec3_sub(point, a));
}

/* Clip a segment against the plane and inward half-planes of a triangle. */
static bool clip_segment_to_triangle(geomlib_vec3 a, geomlib_vec3 b,
                                     const geomlib_vec3 triangle[3], geomlib_vec3 normal,
                                     float tolerance, float *first, float *last) {
    float da = geomlib_vec3_dot(normal, geomlib_vec3_sub(a, triangle[0]));
    float db = geomlib_vec3_dot(normal, geomlib_vec3_sub(b, triangle[0]));
    *first = 0.0f;
    *last = 1.0f;
    if (fabsf(da) > tolerance || fabsf(db) > tolerance) {
        if ((da > tolerance && db > tolerance) || (da < -tolerance && db < -tolerance))
            return false;
        *first = *last = fmaxf(0.0f, fminf(1.0f, da / (da - db)));
    }
    for (int i = 0; i < 3; i++) {
        geomlib_vec3 p = triangle[i], q = triangle[(i + 1) % 3];
        float start = triangle_half_plane(p, q, normal, a);
        float end = triangle_half_plane(p, q, normal, b);
        if (start < -tolerance && end < -tolerance) return false;
        if (start < -tolerance) *first = fmaxf(*first, (start + tolerance) / (start - end));
        if (end < -tolerance) *last = fminf(*last, (start + tolerance) / (start - end));
        if (*first > *last) return false;
    }
    return true;
}

static bool on_shared_vertex(geomlib_vec3 point, const geomlib_vec3 shared[2],
                             int shared_count, float tolerance) {
    if (!shared_count) return false;
    geomlib_vec3 offset = geomlib_vec3_sub(point, shared[0]);
    return geomlib_vec3_length(offset) <= tolerance;
}

static bool triangles_intersect(const geomlib_trimesh *mesh, uint32_t face_a, uint32_t face_b,
                                const geomlib_vec3 a[3], const geomlib_vec3 b[3],
                                float tolerance, geomlib_vec3 *point) {
    const geomlib_vec3 *triangles[2] = {a, b};
    geomlib_vec3 normals[2];
    for (int i = 0; i < 2; i++) {
        const geomlib_vec3 *triangle = triangles[i];
        geomlib_vec3 u = geomlib_vec3_sub(triangle[1], triangle[0]);
        geomlib_vec3 v = geomlib_vec3_sub(triangle[2], triangle[0]);
        geomlib_vec3 normal = geomlib_vec3_cross(u, v);
        float length = geomlib_vec3_length(normal);
        if (!(length > 0.0f)) return false;
        normals[i] = geomlib_vec3_mul(normal, 1.0f / length);
    }

    geomlib_vec3 shared[2];
    int shared_count = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (mesh->indices[face_a * 3 + i] != mesh->indices[face_b * 3 + j]) continue;
            if (shared_count == 2) {
                *point = a[0];  /* Duplicate triangles overlap across their full area. */
                return true;
            }
            shared[shared_count++] = a[i];
            break;
        }
    }

    if (shared_count == 2) {
        geomlib_vec3 edge = geomlib_vec3_sub(shared[1], shared[0]);
        geomlib_vec3 other[2] = {{0, 0, 0}, {0, 0, 0}};
        for (int side = 0; side < 2; side++) {
            for (int i = 0; i < 3; i++) {
                geomlib_vec3 offset = geomlib_vec3_sub(triangles[side][i], shared[0]);
                geomlib_vec3 perpendicular = geomlib_vec3_cross(edge, offset);
                if (geomlib_vec3_length(perpendicular) > geomlib_vec3_length(other[side]))
                    other[side] = perpendicular;
            }
        }
        if (geomlib_vec3_length(geomlib_vec3_cross(normals[0], normals[1])) > 1e-6f ||
            geomlib_vec3_dot(other[0], other[1]) <= 0.0f) return false;
        *point = shared[0];
        return true;
    }

    for (int side = 0; side < 2; side++) {
        const geomlib_vec3 *source = triangles[side], *target = triangles[1 - side];
        for (int i = 0; i < 3; i++) {
            geomlib_vec3 start = source[i], end = source[(i + 1) % 3];
            float first, last;
            if (!clip_segment_to_triangle(start, end, target, normals[1 - side],
                                          tolerance, &first, &last)) continue;
            geomlib_vec3 edge = geomlib_vec3_sub(end, start);
            geomlib_vec3 points[2] = {geomlib_vec3_add(start, geomlib_vec3_mul(edge, first)),
                                      geomlib_vec3_add(start, geomlib_vec3_mul(edge, last))};
            for (int j = 0; j < 2; j++) {
                if (on_shared_vertex(points[j], shared, shared_count, tolerance * 4.0f)) continue;
                *point = points[j];
                return true;
            }
        }
    }
    return false;
}

static bool collect_intersections(const geomlib_bvh *bvh, uint32_t node_index,
                                  uint32_t face, const geomlib_vec3 triangle[3],
                                  geomlib_aabb bounds, geomlib_self_intersection *result,
                                  size_t *capacity) {
    const geomlib_bvh_node *node = &bvh->nodes[node_index];
    if (!geomlib_aabb_intersects(bounds, node->bounds)) return true;
    if (!node->is_leaf) {
        return collect_intersections(bvh, node->data.interior.left_child, face, triangle,
                                     bounds, result, capacity) &&
               collect_intersections(bvh, node->data.interior.right_child, face, triangle,
                                     bounds, result, capacity);
    }

    for (uint32_t i = 0; i < node->data.leaf.primitive_count; i++) {
        uint32_t other = bvh->primitive_indices[node->data.leaf.primitive_offset + i];
        if (other <= face) continue;
        geomlib_vec3 vertices[3];
        geomlib_trimesh_get_triangle_vertices(bvh->mesh, other, &vertices[0], &vertices[1], &vertices[2]);
        geomlib_aabb other_bounds = triangle_bounds(vertices);
        if (!geomlib_aabb_intersects(bounds, other_bounds)) continue;
        float scale = fminf(geomlib_vec3_length(geomlib_vec3_sub(bounds.max, bounds.min)),
                            geomlib_vec3_length(geomlib_vec3_sub(other_bounds.max, other_bounds.min)));
        geomlib_vec3 point;
        if (!triangles_intersect(bvh->mesh, face, other, triangle, vertices, scale * 1e-6f, &point))
            continue;
        if (result->pair_count == *capacity) {
            size_t next = *capacity ? *capacity * 2 : 16;
            if (next <= *capacity || next > SIZE_MAX / sizeof(*result->pairs)) return false;
            geomlib_intersection_pair *pairs = (geomlib_intersection_pair*)realloc(
                result->pairs, next * sizeof(*pairs));
            if (!pairs) return false;
            result->pairs = pairs;
            *capacity = next;
        }
        result->pairs[result->pair_count++] = (geomlib_intersection_pair){face, other, point};
    }
    return true;
}

geomlib_self_intersection geomlib_detect_self_intersections(const geomlib_trimesh *mesh) {
    geomlib_self_intersection result = {0};
    if (!mesh || !mesh->triangle_count || mesh->triangle_count > UINT32_MAX / 3 ||
        !mesh->vertices || !mesh->indices) return result;
    for (size_t i = 0; i < mesh->triangle_count * 3; i++) {
        if (mesh->indices[i] >= mesh->vertex_count) return result;
    }
    geomlib_bvh bvh;
    if (geomlib_bvh_build(&bvh, mesh, GEOMLIB_BVH_QUALITY_MEDIUM) != GEOMLIB_MESH_SUCCESS)
        return result;

    size_t capacity = 0;
    for (uint32_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 triangle[3];
        geomlib_trimesh_get_triangle_vertices(mesh, i, &triangle[0], &triangle[1], &triangle[2]);
        if (!collect_intersections(&bvh, 0, i, triangle, triangle_bounds(triangle), &result, &capacity)) {
            geomlib_self_intersection_free(&result);
            result = (geomlib_self_intersection){0};
            break;
        }
    }
    result.has_intersections = result.pair_count != 0;
    geomlib_bvh_free(&bvh);
    return result;
}

void geomlib_self_intersection_free(geomlib_self_intersection *result) {
    free(result->pairs);
}

/* ============================================================================
 * Face Adjacency
 * ============================================================================ */

geomlib_face_adjacency geomlib_build_face_adjacency(const geomlib_trimesh *mesh) {
    geomlib_face_adjacency result = {0};
    edge_hash_map map = {0};
    if (build_edge_map(mesh, &map, NULL, NULL) && mesh->triangle_count) {
        size_t count = mesh->triangle_count * 3;
        result.neighbors = (uint32_t*)malloc(count * sizeof(uint32_t));
        if (result.neighbors) {
            for (size_t i = 0; i < count; i++) {
                size_t next = i % 3 == 2 ? i - 2 : i + 1;
                edge_key key = make_edge_key(mesh->indices[i], mesh->indices[next]);
                edge_hash_entry *edge = edge_hash_map_find(&map, key);
                result.neighbors[i] = edge->face_count == 2
                    ? edge->faces[edge->faces[0] == i / 3] : UINT32_MAX;
            }
        }
    }
    edge_hash_map_free(&map);
    return result;
}

void geomlib_face_adjacency_free(geomlib_face_adjacency *adjacency) {
    free(adjacency->neighbors);
}

/* ============================================================================
 * Topological Operations
 * ============================================================================ */

int geomlib_remove_isolated_vertices(geomlib_trimesh *mesh) {
    if (!mesh || mesh->vertex_count == 0) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Mark vertices that are referenced by faces */
    bool *used = (bool*)calloc(mesh->vertex_count, sizeof(bool));
    if (!used) {
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (size_t i = 0; i < mesh->triangle_count * 3; i++) {
        if (mesh->indices[i] < mesh->vertex_count) {
            used[mesh->indices[i]] = true;
        }
    }

    /* Count isolated vertices */
    size_t isolated_count = 0;
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        if (!used[i]) {
            isolated_count++;
        }
    }

    if (isolated_count == 0) {
        free(used);
        return GEOMLIB_MESH_SUCCESS;
    }

    /* Build vertex remap table */
    uint32_t *remap = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    if (!remap) {
        free(used);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    size_t new_index = 0;
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        if (used[i]) {
            remap[i] = (uint32_t)new_index;
            new_index++;
        } else {
            remap[i] = UINT32_MAX;
        }
    }

    /* Create new vertex array */
    size_t new_vertex_count = mesh->vertex_count - isolated_count;
    geomlib_vec3 *new_vertices = (geomlib_vec3*)malloc(new_vertex_count * sizeof(geomlib_vec3));
    if (!new_vertices) {
        free(used);
        free(remap);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (size_t i = 0; i < mesh->vertex_count; i++) {
        if (used[i]) {
            new_vertices[remap[i]] = mesh->vertices[i];
        }
    }

    /* Update indices */
    for (size_t i = 0; i < mesh->triangle_count * 3; i++) {
        mesh->indices[i] = remap[mesh->indices[i]];
    }

    /* Replace vertices */
    free(mesh->vertices);
    mesh->vertices = new_vertices;
    mesh->vertex_count = new_vertex_count;
    mesh->vertex_capacity = new_vertex_count;

    free(used);
    free(remap);

    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_remove_duplicate_faces(geomlib_trimesh *mesh) {
    if (!mesh || mesh->triangle_count == 0) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Mark duplicate faces */
    bool *keep = (bool*)malloc(mesh->triangle_count * sizeof(bool));
    if (!keep) {
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        keep[i] = true;
    }

    /* Check each pair (O(n²) but simple and correct) */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        if (!keep[i]) continue;

        uint32_t v0 = mesh->indices[i * 3 + 0];
        uint32_t v1 = mesh->indices[i * 3 + 1];
        uint32_t v2 = mesh->indices[i * 3 + 2];

        for (size_t j = i + 1; j < mesh->triangle_count; j++) {
            if (!keep[j]) continue;

            uint32_t u0 = mesh->indices[j * 3 + 0];
            uint32_t u1 = mesh->indices[j * 3 + 1];
            uint32_t u2 = mesh->indices[j * 3 + 2];

            /* Check if same triangle (any rotation/permutation) */
            bool same = false;
            if ((v0 == u0 && v1 == u1 && v2 == u2) ||
                (v0 == u1 && v1 == u2 && v2 == u0) ||
                (v0 == u2 && v1 == u0 && v2 == u1)) {
                same = true;
            }

            if (same) {
                keep[j] = false;
            }
        }
    }

    /* Count faces to keep */
    size_t new_count = 0;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        if (keep[i]) {
            new_count++;
        }
    }

    if (new_count == mesh->triangle_count) {
        free(keep);
        return GEOMLIB_MESH_SUCCESS;
    }

    /* Compact indices */
    size_t write_idx = 0;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        if (keep[i]) {
            if (write_idx != i) {
                mesh->indices[write_idx * 3 + 0] = mesh->indices[i * 3 + 0];
                mesh->indices[write_idx * 3 + 1] = mesh->indices[i * 3 + 1];
                mesh->indices[write_idx * 3 + 2] = mesh->indices[i * 3 + 2];
            }
            write_idx++;
        }
    }

    mesh->triangle_count = new_count;

    free(keep);
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_split_nonmanifold_vertices(geomlib_trimesh *mesh) {
    if (!mesh || !mesh->vertex_count || !mesh->triangle_count) return 0;
    edge_hash_map map = {0};
    uint32_t *fans = build_corner_fans(mesh, &map);
    edge_hash_map_free(&map);
    if (!fans) return -GEOMLIB_MESH_ERROR_MEMORY;

    size_t original_count = mesh->vertex_count;
    size_t corner_count = mesh->triangle_count * 3;
    uint32_t *remap = (uint32_t*)malloc(corner_count * sizeof(uint32_t));
    bool *has_fan = (bool*)calloc(original_count, sizeof(bool));
    int result = -GEOMLIB_MESH_ERROR_MEMORY;
    if (!remap || !has_fan) goto finish;

    size_t new_count = original_count;
    for (size_t i = 0; i < corner_count; i++) {
        if (fans[i] != i) continue;
        uint32_t vertex = mesh->indices[i];
        if (!has_fan[vertex]) {
            remap[i] = vertex;
            has_fan[vertex] = true;
        } else {
            if (new_count >= UINT32_MAX || new_count - original_count >= INT_MAX) goto finish;
            remap[i] = (uint32_t)new_count++;
        }
    }
    if (new_count > SIZE_MAX / sizeof(geomlib_vec3)) goto finish;
    if (geomlib_trimesh_reserve_vertices(mesh, new_count) != GEOMLIB_MESH_SUCCESS) goto finish;
    for (size_t i = 0; i < corner_count; i++) {
        if (fans[i] != i || remap[i] < original_count) continue;
        uint32_t from = mesh->indices[i], to = remap[i];
        mesh->vertices[to] = mesh->vertices[from];
        if (mesh->normals) mesh->normals[to] = mesh->normals[from];
        if (mesh->texcoords) mesh->texcoords[to] = mesh->texcoords[from];
    }
    for (size_t i = 0; i < corner_count; i++) mesh->indices[i] = remap[fans[i]];
    mesh->vertex_count = new_count;
    result = (int)(new_count - original_count);

finish:
    free(fans);
    free(remap);
    free(has_fan);
    return result;
}

bool geomlib_is_consistently_oriented(const geomlib_trimesh *mesh) {
    if (!mesh || !mesh->triangle_count) return true;
    edge_hash_map map = {0};
    bool consistent = build_edge_map(mesh, &map, NULL, NULL);
    for (size_t i = 0; consistent && i < map.bucket_count; i++) {
        for (edge_hash_entry *edge = map.buckets[i]; edge; edge = edge->next) {
            if (edge->face_count > 2 || (edge->face_count == 2 && edge->forward_count != 1)) {
                consistent = false;
                break;
            }
        }
    }
    edge_hash_map_free(&map);
    return consistent;
}

int geomlib_orient_consistently(geomlib_trimesh *mesh) {
    if (!mesh || !mesh->triangle_count) return GEOMLIB_MESH_ERROR_INVALID;
    edge_hash_map map = {0};
    uint8_t *orientation = (uint8_t*)calloc(mesh->triangle_count, sizeof(uint8_t));
    uint32_t *queue = (uint32_t*)malloc(mesh->triangle_count * sizeof(uint32_t));
    int result = GEOMLIB_MESH_ERROR_MEMORY;
    if (!orientation || !queue || !build_edge_map(mesh, &map, NULL, NULL)) goto finish;

    result = GEOMLIB_MESH_ERROR_INVALID;
    for (size_t seed = 0; seed < mesh->triangle_count; seed++) {
        if (orientation[seed]) continue;
        size_t front = 0, end = 0;
        queue[end++] = (uint32_t)seed;
        orientation[seed] = 1;
        while (front < end) {
            uint32_t face = queue[front++];
            for (int i = 0; i < 3; i++) {
                uint32_t a = mesh->indices[face * 3 + i];
                uint32_t b = mesh->indices[face * 3 + (i + 1) % 3];
                edge_hash_entry *edge = edge_hash_map_find(&map, make_edge_key(a, b));
                if (edge->face_count == 1) continue;
                if (edge->face_count > 2) goto finish;
                uint32_t neighbor = edge->faces[edge->faces[0] == face];
                uint8_t required = edge->forward_count == 1 ? orientation[face] : 3 - orientation[face];
                if (!orientation[neighbor]) {
                    orientation[neighbor] = required;
                    queue[end++] = neighbor;
                } else if (orientation[neighbor] != required) {
                    goto finish;
                }
            }
        }
    }

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        if (orientation[i] != 2) continue;
        uint32_t temp = mesh->indices[i * 3 + 1];
        mesh->indices[i * 3 + 1] = mesh->indices[i * 3 + 2];
        mesh->indices[i * 3 + 2] = temp;
    }
    result = GEOMLIB_MESH_SUCCESS;

finish:
    free(orientation);
    free(queue);
    edge_hash_map_free(&map);
    return result;
}

int geomlib_reverse_orientation(geomlib_trimesh *mesh) {
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t temp = mesh->indices[i * 3 + 1];
        mesh->indices[i * 3 + 1] = mesh->indices[i * 3 + 2];
        mesh->indices[i * 3 + 2] = temp;
    }
    return GEOMLIB_MESH_SUCCESS;
}
