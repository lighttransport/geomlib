/*
 * geomlib - Mesh Topology Analysis Implementation
 */

#include "geomlib/geomlib_topology.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
    uint32_t faces[4];  /* Up to 4 faces (non-manifold if >2) */
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
    for (size_t i = 0; i < map->bucket_count; i++) {
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

geomlib_edge_mesh geomlib_build_edge_mesh(const geomlib_trimesh *mesh) {
    geomlib_edge_mesh edge_mesh;

    /* Build edge hash map */
    size_t estimated_edges = mesh->triangle_count * 3 / 2;  /* Euler: E ≈ 1.5F for closed */
    edge_hash_map map = edge_hash_map_create(estimated_edges * 2);

    /* Insert all edges */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t v0 = mesh->indices[i * 3 + 0];
        uint32_t v1 = mesh->indices[i * 3 + 1];
        uint32_t v2 = mesh->indices[i * 3 + 2];

        edge_key edges[3] = {
            make_edge_key(v0, v1),
            make_edge_key(v1, v2),
            make_edge_key(v2, v0)
        };

        for (int j = 0; j < 3; j++) {
            edge_hash_entry *entry = edge_hash_map_insert(&map, edges[j]);
            if (entry->face_count < 4) {
                entry->faces[entry->face_count++] = (uint32_t)i;
            }
        }
    }

    /* Build edge array */
    edge_mesh.edge_count = map.entry_count;
    edge_mesh.edges = (geomlib_edge*)malloc(edge_mesh.edge_count * sizeof(geomlib_edge));

    size_t edge_idx = 0;
    for (size_t i = 0; i < map.bucket_count; i++) {
        edge_hash_entry *entry = map.buckets[i];
        while (entry) {
            edge_mesh.edges[edge_idx].v0 = entry->key.v0;
            edge_mesh.edges[edge_idx].v1 = entry->key.v1;
            edge_mesh.edges[edge_idx].face0 = entry->faces[0];
            edge_mesh.edges[edge_idx].face1 = entry->face_count > 1 ? entry->faces[1] : UINT32_MAX;
            edge_idx++;
            entry = entry->next;
        }
    }

    /* Build vertex edge lists */
    edge_mesh.vertex_edges = (uint32_t**)calloc(mesh->vertex_count, sizeof(uint32_t*));
    edge_mesh.vertex_edge_count = (uint32_t*)calloc(mesh->vertex_count, sizeof(uint32_t));

    /* Count edges per vertex */
    for (size_t i = 0; i < edge_mesh.edge_count; i++) {
        edge_mesh.vertex_edge_count[edge_mesh.edges[i].v0]++;
        edge_mesh.vertex_edge_count[edge_mesh.edges[i].v1]++;
    }

    /* Allocate and fill */
    uint32_t *counts = (uint32_t*)calloc(mesh->vertex_count, sizeof(uint32_t));
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        edge_mesh.vertex_edges[i] = (uint32_t*)malloc(
            edge_mesh.vertex_edge_count[i] * sizeof(uint32_t));
    }

    for (size_t i = 0; i < edge_mesh.edge_count; i++) {
        uint32_t v0 = edge_mesh.edges[i].v0;
        uint32_t v1 = edge_mesh.edges[i].v1;
        edge_mesh.vertex_edges[v0][counts[v0]++] = (uint32_t)i;
        edge_mesh.vertex_edges[v1][counts[v1]++] = (uint32_t)i;
    }

    free(counts);
    edge_hash_map_free(&map);

    return edge_mesh;
}

void geomlib_edge_mesh_free(geomlib_edge_mesh *edge_mesh) {
    free(edge_mesh->edges);
    for (size_t i = 0; i < edge_mesh->edge_count; i++) {
        free(edge_mesh->vertex_edges[i]);
    }
    free(edge_mesh->vertex_edges);
    free(edge_mesh->vertex_edge_count);
}

/* ============================================================================
 * Euler Characteristic
 * ============================================================================ */

geomlib_euler_characteristic geomlib_compute_euler_characteristic(
    const geomlib_trimesh *mesh
) {
    geomlib_euler_characteristic result;

    result.vertices = (int32_t)mesh->vertex_count;
    result.faces = (int32_t)mesh->triangle_count;

    /* Build edge mesh to count edges */
    geomlib_edge_mesh edge_mesh = geomlib_build_edge_mesh(mesh);
    result.edges = (int32_t)edge_mesh.edge_count;

    /* Count boundary edges */
    int32_t boundary_edges = 0;
    for (size_t i = 0; i < edge_mesh.edge_count; i++) {
        if (edge_mesh.edges[i].face1 == UINT32_MAX) {
            boundary_edges++;
        }
    }

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

geomlib_boundary_loops geomlib_extract_boundary_loops(
    const geomlib_trimesh *mesh
) {
    geomlib_boundary_loops result;
    result.loops = NULL;
    result.loop_count = 0;

    geomlib_edge_mesh edge_mesh = geomlib_build_edge_mesh(mesh);

    /* Mark boundary edges */
    bool *visited = (bool*)calloc(edge_mesh.edge_count, sizeof(bool));

    /* Dynamic array for loops */
    size_t loop_capacity = 4;
    result.loops = (geomlib_boundary_loop*)malloc(loop_capacity * sizeof(geomlib_boundary_loop));

    for (size_t start_edge = 0; start_edge < edge_mesh.edge_count; start_edge++) {
        if (visited[start_edge] || edge_mesh.edges[start_edge].face1 != UINT32_MAX) {
            continue;  /* Not a boundary edge or already visited */
        }

        /* Start new loop */
        size_t vertex_capacity = 16;
        uint32_t *loop_vertices = (uint32_t*)malloc(vertex_capacity * sizeof(uint32_t));
        size_t vertex_count = 0;
        float total_length = 0.0f;

        /* Walk the boundary */
        uint32_t current_edge = (uint32_t)start_edge;
        uint32_t current_vertex = edge_mesh.edges[current_edge].v0;

        do {
            /* Add vertex */
            if (vertex_count >= vertex_capacity) {
                vertex_capacity *= 2;
                loop_vertices = (uint32_t*)realloc(loop_vertices, vertex_capacity * sizeof(uint32_t));
            }
            loop_vertices[vertex_count++] = current_vertex;

            /* Mark edge as visited */
            visited[current_edge] = true;

            /* Find next boundary edge */
            uint32_t next_vertex = (edge_mesh.edges[current_edge].v0 == current_vertex) ?
                edge_mesh.edges[current_edge].v1 : edge_mesh.edges[current_edge].v0;

            /* Add edge length */
            geomlib_vec3 v0 = mesh->vertices[current_vertex];
            geomlib_vec3 v1 = mesh->vertices[next_vertex];
            total_length += geomlib_vec3_length(geomlib_vec3_sub(v1, v0));

            /* Find next edge */
            bool found = false;
            for (uint32_t i = 0; i < edge_mesh.vertex_edge_count[next_vertex]; i++) {
                uint32_t candidate = edge_mesh.vertex_edges[next_vertex][i];
                if (candidate != current_edge &&
                    edge_mesh.edges[candidate].face1 == UINT32_MAX &&
                    !visited[candidate]) {
                    current_edge = candidate;
                    current_vertex = next_vertex;
                    found = true;
                    break;
                }
            }

            if (!found) {
                break;  /* Loop complete or error */
            }
        } while (current_edge != start_edge && vertex_count < mesh->vertex_count);

        /* Store loop */
        if (result.loop_count >= loop_capacity) {
            loop_capacity *= 2;
            result.loops = (geomlib_boundary_loop*)realloc(
                result.loops, loop_capacity * sizeof(geomlib_boundary_loop));
        }

        result.loops[result.loop_count].vertices = loop_vertices;
        result.loops[result.loop_count].vertex_count = vertex_count;
        result.loops[result.loop_count].length = total_length;
        result.loop_count++;
    }

    free(visited);
    geomlib_edge_mesh_free(&edge_mesh);

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

geomlib_manifold_analysis geomlib_analyze_manifold(
    const geomlib_trimesh *mesh
) {
    geomlib_manifold_analysis analysis;
    analysis.is_edge_manifold = true;
    analysis.is_vertex_manifold = true;
    analysis.is_oriented = true;
    analysis.nonmanifold_elements = NULL;
    analysis.nonmanifold_count = 0;

    geomlib_edge_mesh edge_mesh = geomlib_build_edge_mesh(mesh);

    /* Check edge manifold (each edge has at most 2 faces) */
    size_t nm_capacity = 16;
    geomlib_nonmanifold_element *nm_elements = (geomlib_nonmanifold_element*)malloc(
        nm_capacity * sizeof(geomlib_nonmanifold_element));

    size_t boundary_count = 0;
    for (size_t i = 0; i < edge_mesh.edge_count; i++) {
        uint32_t face_count = (edge_mesh.edges[i].face1 == UINT32_MAX) ? 1 : 2;

        /* Count actual faces (check for >2) */
        edge_key key = make_edge_key(edge_mesh.edges[i].v0, edge_mesh.edges[i].v1);
        size_t actual_count = 0;
        for (size_t j = 0; j < mesh->triangle_count; j++) {
            for (int k = 0; k < 3; k++) {
                uint32_t v0 = mesh->indices[j * 3 + k];
                uint32_t v1 = mesh->indices[j * 3 + ((k + 1) % 3)];
                edge_key test = make_edge_key(v0, v1);
                if (test.v0 == key.v0 && test.v1 == key.v1) {
                    actual_count++;
                }
            }
        }

        if (actual_count > 2) {
            analysis.is_edge_manifold = false;

            if (analysis.nonmanifold_count >= nm_capacity) {
                nm_capacity *= 2;
                nm_elements = (geomlib_nonmanifold_element*)realloc(
                    nm_elements, nm_capacity * sizeof(geomlib_nonmanifold_element));
            }

            nm_elements[analysis.nonmanifold_count].type = GEOMLIB_NONMANIFOLD_EDGE;
            nm_elements[analysis.nonmanifold_count].index = (uint32_t)i;
            nm_elements[analysis.nonmanifold_count].face_count = (uint32_t)actual_count;
            analysis.nonmanifold_count++;
        }

        if (face_count == 1) {
            boundary_count++;
        }
    }

    analysis.is_closed = (boundary_count == 0);
    analysis.is_manifold = analysis.is_edge_manifold && analysis.is_vertex_manifold;
    analysis.nonmanifold_elements = nm_elements;

    geomlib_edge_mesh_free(&edge_mesh);

    return analysis;
}

void geomlib_manifold_analysis_free(geomlib_manifold_analysis *analysis) {
    free(analysis->nonmanifold_elements);
}

/* ============================================================================
 * Vertex Adjacency
 * ============================================================================ */

geomlib_vertex_adjacency geomlib_build_vertex_adjacency(
    const geomlib_trimesh *mesh
) {
    geomlib_vertex_adjacency adjacency;

    adjacency.neighbors = (uint32_t**)calloc(mesh->vertex_count, sizeof(uint32_t*));
    adjacency.neighbor_counts = (uint32_t*)calloc(mesh->vertex_count, sizeof(uint32_t));
    adjacency.faces = (uint32_t**)calloc(mesh->vertex_count, sizeof(uint32_t*));
    adjacency.face_counts = (uint32_t*)calloc(mesh->vertex_count, sizeof(uint32_t));

    /* Count adjacencies */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        for (int j = 0; j < 3; j++) {
            uint32_t v = mesh->indices[i * 3 + j];
            adjacency.neighbor_counts[v] += 2;  /* Two neighbors per face */
            adjacency.face_counts[v]++;
        }
    }

    /* Allocate */
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        adjacency.neighbors[i] = (uint32_t*)malloc(adjacency.neighbor_counts[i] * sizeof(uint32_t));
        adjacency.faces[i] = (uint32_t*)malloc(adjacency.face_counts[i] * sizeof(uint32_t));
    }

    /* Fill */
    uint32_t *n_counts = (uint32_t*)calloc(mesh->vertex_count, sizeof(uint32_t));
    uint32_t *f_counts = (uint32_t*)calloc(mesh->vertex_count, sizeof(uint32_t));

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t v0 = mesh->indices[i * 3 + 0];
        uint32_t v1 = mesh->indices[i * 3 + 1];
        uint32_t v2 = mesh->indices[i * 3 + 2];

        adjacency.neighbors[v0][n_counts[v0]++] = v1;
        adjacency.neighbors[v0][n_counts[v0]++] = v2;
        adjacency.neighbors[v1][n_counts[v1]++] = v0;
        adjacency.neighbors[v1][n_counts[v1]++] = v2;
        adjacency.neighbors[v2][n_counts[v2]++] = v0;
        adjacency.neighbors[v2][n_counts[v2]++] = v1;

        adjacency.faces[v0][f_counts[v0]++] = (uint32_t)i;
        adjacency.faces[v1][f_counts[v1]++] = (uint32_t)i;
        adjacency.faces[v2][f_counts[v2]++] = (uint32_t)i;
    }

    free(n_counts);
    free(f_counts);

    return adjacency;
}

void geomlib_vertex_adjacency_free(geomlib_vertex_adjacency *adjacency) {
    for (size_t i = 0; i < adjacency->neighbor_counts[0]; i++) {  /* Assumes all allocated */
        free(adjacency->neighbors[i]);
        free(adjacency->faces[i]);
    }
    free(adjacency->neighbors);
    free(adjacency->neighbor_counts);
    free(adjacency->faces);
    free(adjacency->face_counts);
}

/* ============================================================================
 * Valence Statistics
 * ============================================================================ */

geomlib_valence_stats geomlib_compute_valence_stats(
    const geomlib_trimesh *mesh
) {
    geomlib_valence_stats stats;

    geomlib_vertex_adjacency adjacency = geomlib_build_vertex_adjacency(mesh);

    stats.valences = (uint32_t*)malloc(mesh->vertex_count * sizeof(uint32_t));
    stats.min_valence = UINT32_MAX;
    stats.max_valence = 0;
    uint32_t total_valence = 0;

    /* Remove duplicates and count unique neighbors */
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        /* Use simple O(n²) duplicate removal for small neighborhoods */
        uint32_t unique_count = 0;
        for (uint32_t j = 0; j < adjacency.neighbor_counts[i]; j++) {
            bool is_duplicate = false;
            uint32_t neighbor = adjacency.neighbors[i][j];
            for (uint32_t k = 0; k < j; k++) {
                if (adjacency.neighbors[i][k] == neighbor) {
                    is_duplicate = true;
                    break;
                }
            }
            if (!is_duplicate) {
                unique_count++;
            }
        }

        stats.valences[i] = unique_count;
        if (unique_count < stats.min_valence) stats.min_valence = unique_count;
        if (unique_count > stats.max_valence) stats.max_valence = unique_count;
        total_valence += unique_count;
    }

    stats.avg_valence = (float)total_valence / (float)mesh->vertex_count;

    /* Build histogram */
    stats.histogram_size = stats.max_valence + 1;
    stats.histogram = (uint32_t*)calloc(stats.histogram_size, sizeof(uint32_t));
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        stats.histogram[stats.valences[i]]++;
    }

    geomlib_vertex_adjacency_free(&adjacency);

    return stats;
}

void geomlib_valence_stats_free(geomlib_valence_stats *stats) {
    free(stats->valences);
    free(stats->histogram);
}

/* ============================================================================
 * Self-Intersection Detection (Basic BVH-based)
 * ============================================================================ */

geomlib_self_intersection geomlib_detect_self_intersections(
    const geomlib_trimesh *mesh
) {
    /* TODO: Implement using BVH for efficiency */
    geomlib_self_intersection result;
    result.has_intersections = false;
    result.pairs = NULL;
    result.pair_count = 0;
    return result;
}

void geomlib_self_intersection_free(geomlib_self_intersection *result) {
    free(result->pairs);
}

/* ============================================================================
 * Face Adjacency
 * ============================================================================ */

geomlib_face_adjacency geomlib_build_face_adjacency(
    const geomlib_trimesh *mesh
) {
    geomlib_face_adjacency adjacency;
    adjacency.neighbors = (uint32_t*)malloc(mesh->triangle_count * 3 * sizeof(uint32_t));

    /* Initialize to boundary */
    for (size_t i = 0; i < mesh->triangle_count * 3; i++) {
        adjacency.neighbors[i] = UINT32_MAX;
    }

    /* Build edge map */
    edge_hash_map map = edge_hash_map_create(mesh->triangle_count * 3);

    for (size_t i = 0; i < mesh->triangle_count; i++) {
        for (int j = 0; j < 3; j++) {
            uint32_t v0 = mesh->indices[i * 3 + j];
            uint32_t v1 = mesh->indices[i * 3 + ((j + 1) % 3)];
            edge_key key = make_edge_key(v0, v1);

            edge_hash_entry *entry = edge_hash_map_insert(&map, key);
            if (entry->face_count < 2) {
                entry->faces[entry->face_count++] = (uint32_t)i;
            }
        }
    }

    /* Fill adjacency */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        for (int j = 0; j < 3; j++) {
            uint32_t v0 = mesh->indices[i * 3 + j];
            uint32_t v1 = mesh->indices[i * 3 + ((j + 1) % 3)];
            edge_key key = make_edge_key(v0, v1);

            edge_hash_entry *entry = edge_hash_map_find(&map, key);
            if (entry && entry->face_count == 2) {
                uint32_t other = (entry->faces[0] == i) ? entry->faces[1] : entry->faces[0];
                adjacency.neighbors[i * 3 + j] = other;
            }
        }
    }

    edge_hash_map_free(&map);

    return adjacency;
}

void geomlib_face_adjacency_free(geomlib_face_adjacency *adjacency) {
    free(adjacency->neighbors);
}

/* ============================================================================
 * Topological Operations (Stubs)
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
    /* Complex operation - stub for now */
    /* TODO: Implement proper nonmanifold vertex splitting */
    (void)mesh;
    return GEOMLIB_MESH_SUCCESS;
}

bool geomlib_is_consistently_oriented(const geomlib_trimesh *mesh) {
    if (!mesh || mesh->triangle_count == 0) {
        return true;
    }

    /* Build edge mesh to check orientation consistency */
    geomlib_edge_mesh edges = geomlib_build_edge_mesh(mesh);

    bool is_consistent = true;

    /* Check each edge - if it has 2 adjacent faces, they should have opposite orientations */
    for (size_t i = 0; i < edges.edge_count; i++) {
        if (edges.edges[i].face1 == UINT32_MAX) {
            continue;  /* Boundary edge */
        }

        uint32_t v0 = edges.edges[i].v0;
        uint32_t v1 = edges.edges[i].v1;
        uint32_t f0 = edges.edges[i].face0;
        uint32_t f1 = edges.edges[i].face1;

        /* Find edge in face0 */
        bool f0_forward = false;
        for (int j = 0; j < 3; j++) {
            uint32_t a = mesh->indices[f0 * 3 + j];
            uint32_t b = mesh->indices[f0 * 3 + ((j + 1) % 3)];
            if ((a == v0 && b == v1)) {
                f0_forward = true;
                break;
            }
        }

        /* Find edge in face1 */
        bool f1_forward = false;
        for (int j = 0; j < 3; j++) {
            uint32_t a = mesh->indices[f1 * 3 + j];
            uint32_t b = mesh->indices[f1 * 3 + ((j + 1) % 3)];
            if ((a == v0 && b == v1)) {
                f1_forward = true;
                break;
            }
        }

        /* Edges should have opposite orientations */
        if (f0_forward == f1_forward) {
            is_consistent = false;
            break;
        }
    }

    geomlib_edge_mesh_free(&edges);

    return is_consistent;
}

int geomlib_orient_consistently(geomlib_trimesh *mesh) {
    if (!mesh || mesh->triangle_count == 0) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Build face adjacency */
    geomlib_edge_mesh edges = geomlib_build_edge_mesh(mesh);

    /* Track which faces have been oriented */
    bool *oriented = (bool*)calloc(mesh->triangle_count, sizeof(bool));
    bool *flipped = (bool*)calloc(mesh->triangle_count, sizeof(bool));

    if (!oriented || !flipped) {
        free(oriented);
        free(flipped);
        geomlib_edge_mesh_free(&edges);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    /* BFS to propagate orientation */
    uint32_t *queue = (uint32_t*)malloc(mesh->triangle_count * sizeof(uint32_t));
    if (!queue) {
        free(oriented);
        free(flipped);
        geomlib_edge_mesh_free(&edges);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    /* Process each connected component */
    for (size_t seed = 0; seed < mesh->triangle_count; seed++) {
        if (oriented[seed]) continue;

        /* Start BFS from this seed */
        size_t queue_head = 0;
        size_t queue_tail = 0;

        queue[queue_tail++] = (uint32_t)seed;
        oriented[seed] = true;
        flipped[seed] = false;

        while (queue_head < queue_tail) {
            uint32_t current = queue[queue_head++];

            /* Find adjacent faces via edges */
            for (size_t i = 0; i < edges.edge_count; i++) {
                uint32_t f0 = edges.edges[i].face0;
                uint32_t f1 = edges.edges[i].face1;

                if (f1 == UINT32_MAX) continue;  /* Boundary */

                uint32_t neighbor = UINT32_MAX;
                if (f0 == current && !oriented[f1]) {
                    neighbor = f1;
                } else if (f1 == current && !oriented[f0]) {
                    neighbor = f0;
                }

                if (neighbor == UINT32_MAX) continue;

                /* Check edge orientation between current and neighbor */
                uint32_t v0 = edges.edges[i].v0;
                uint32_t v1 = edges.edges[i].v1;

                /* Find edge direction in current face */
                bool current_forward = false;
                for (int j = 0; j < 3; j++) {
                    uint32_t a = mesh->indices[current * 3 + j];
                    uint32_t b = mesh->indices[current * 3 + ((j + 1) % 3)];
                    if (a == v0 && b == v1) {
                        current_forward = true;
                        break;
                    } else if (a == v1 && b == v0) {
                        current_forward = false;
                        break;
                    }
                }

                /* Find edge direction in neighbor face */
                bool neighbor_forward = false;
                for (int j = 0; j < 3; j++) {
                    uint32_t a = mesh->indices[neighbor * 3 + j];
                    uint32_t b = mesh->indices[neighbor * 3 + ((j + 1) % 3)];
                    if (a == v0 && b == v1) {
                        neighbor_forward = true;
                        break;
                    } else if (a == v1 && b == v0) {
                        neighbor_forward = false;
                        break;
                    }
                }

                /* If orientations match, neighbor needs to be flipped */
                oriented[neighbor] = true;
                if (current_forward == neighbor_forward) {
                    flipped[neighbor] = !flipped[current];
                } else {
                    flipped[neighbor] = flipped[current];
                }

                queue[queue_tail++] = neighbor;
            }
        }
    }

    /* Apply flips */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        if (flipped[i]) {
            uint32_t temp = mesh->indices[i * 3 + 1];
            mesh->indices[i * 3 + 1] = mesh->indices[i * 3 + 2];
            mesh->indices[i * 3 + 2] = temp;
        }
    }

    free(oriented);
    free(flipped);
    free(queue);
    geomlib_edge_mesh_free(&edges);

    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_reverse_orientation(geomlib_trimesh *mesh) {
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t temp = mesh->indices[i * 3 + 1];
        mesh->indices[i * 3 + 1] = mesh->indices[i * 3 + 2];
        mesh->indices[i * 3 + 2] = temp;
    }
    return GEOMLIB_MESH_SUCCESS;
}
