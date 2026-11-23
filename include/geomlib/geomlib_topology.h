/*
 * geomlib - Mesh Topology Analysis
 * Genus computation, connectivity analysis, manifold detection
 */

#ifndef GEOMLIB_TOPOLOGY_H
#define GEOMLIB_TOPOLOGY_H

#include "geomlib_linalg.h"
#include "geomlib_mesh.h"
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Euler Characteristic & Genus
 * ============================================================================ */

/* Topological characteristics */
typedef struct geomlib_euler_characteristic {
    int32_t vertices;           /* V - number of vertices */
    int32_t edges;              /* E - number of edges */
    int32_t faces;              /* F - number of faces */
    int32_t chi;                /* χ = V - E + F (Euler characteristic) */
    int32_t genus;              /* g = (2 - χ - b) / 2 for orientable surfaces */
    int32_t boundaries;         /* b - number of boundary loops */
    int32_t components;         /* Number of connected components */
} geomlib_euler_characteristic;

/* Compute Euler characteristic and genus */
geomlib_euler_characteristic geomlib_compute_euler_characteristic(
    const geomlib_trimesh *mesh
);

/* ============================================================================
 * Edge Representation
 * ============================================================================ */

/* Topology half-edge structure */
typedef struct geomlib_topo_halfedge {
    uint32_t vertex;            /* Origin vertex */
    uint32_t face;              /* Adjacent face */
    uint32_t next;              /* Next halfedge in face */
    uint32_t prev;              /* Previous halfedge in face */
    uint32_t twin;              /* Opposite halfedge (UINT32_MAX if boundary) */
} geomlib_topo_halfedge;

/* Edge (undirected) */
typedef struct geomlib_edge {
    uint32_t v0, v1;            /* Vertex indices (v0 < v1) */
    uint32_t face0, face1;      /* Adjacent faces (face1 = UINT32_MAX if boundary) */
} geomlib_edge;

/* Edge mesh representation */
typedef struct geomlib_edge_mesh {
    geomlib_edge *edges;
    size_t edge_count;
    size_t vertex_count;

    /* Per-vertex edge lists */
    uint32_t **vertex_edges;   /* Edges adjacent to each vertex */
    uint32_t *vertex_edge_count;
} geomlib_edge_mesh;

/* Build edge mesh from triangle mesh */
geomlib_edge_mesh geomlib_build_edge_mesh(const geomlib_trimesh *mesh);

/* Free edge mesh */
void geomlib_edge_mesh_free(geomlib_edge_mesh *edge_mesh);

/* ============================================================================
 * Boundary Detection
 * ============================================================================ */

/* Boundary loop */
typedef struct geomlib_boundary_loop {
    uint32_t *vertices;         /* Ordered vertices around boundary */
    size_t vertex_count;
    float length;               /* Total length of boundary */
} geomlib_boundary_loop;

/* Boundary loops result */
typedef struct geomlib_boundary_loops {
    geomlib_boundary_loop *loops;
    size_t loop_count;
} geomlib_boundary_loops;

/* Extract closed boundary cycles whose vertices each meet two boundary edges. */
geomlib_boundary_loops geomlib_extract_boundary_loops(
    const geomlib_trimesh *mesh
);

/* Free boundary loops */
void geomlib_boundary_loops_free(geomlib_boundary_loops *loops);

/* Check if vertex is on boundary */
bool geomlib_is_boundary_vertex(
    const geomlib_trimesh *mesh,
    uint32_t vertex_idx
);

/* ============================================================================
 * Connected Components
 * ============================================================================ */

/* Connected component info */
typedef struct geomlib_component_info {
    uint32_t *vertex_labels;    /* Component ID for each vertex */
    uint32_t *face_labels;      /* Component ID for each face */
    uint32_t component_count;   /* Number of components */

    /* Per-component statistics */
    uint32_t *component_vertex_counts;
    uint32_t *component_face_counts;
} geomlib_component_info;

/* Find connected components */
geomlib_component_info geomlib_find_connected_components(
    const geomlib_trimesh *mesh
);

/* Free component info */
void geomlib_component_info_free(geomlib_component_info *info);

/* Extract largest component */
int geomlib_extract_largest_component(
    const geomlib_trimesh *mesh,
    geomlib_trimesh *output
);

/* ============================================================================
 * Manifold Detection
 * ============================================================================ */

/* Non-manifold element types */
typedef enum geomlib_nonmanifold_type {
    GEOMLIB_NONMANIFOLD_EDGE,   /* Edge shared by >2 faces */
    GEOMLIB_NONMANIFOLD_VERTEX  /* Vertex with non-disk neighborhood */
} geomlib_nonmanifold_type;

/* Non-manifold element */
typedef struct geomlib_nonmanifold_element {
    geomlib_nonmanifold_type type;
    uint32_t index;             /* Vertex or edge index */
    uint32_t face_count;        /* For edges: number of adjacent faces */
} geomlib_nonmanifold_element;

/* Non-manifold detection result */
typedef struct geomlib_manifold_analysis {
    bool is_manifold;
    bool is_edge_manifold;
    bool is_vertex_manifold;
    bool is_oriented;
    bool is_closed;             /* No boundary */

    geomlib_nonmanifold_element *nonmanifold_elements;
    size_t nonmanifold_count;
} geomlib_manifold_analysis;

/* Check if mesh is manifold */
geomlib_manifold_analysis geomlib_analyze_manifold(
    const geomlib_trimesh *mesh
);

/* Free manifold analysis */
void geomlib_manifold_analysis_free(geomlib_manifold_analysis *analysis);

/* ============================================================================
 * Self-Intersection Detection
 * ============================================================================ */

/* Intersection pair */
typedef struct geomlib_intersection_pair {
    uint32_t face0;             /* First intersecting face */
    uint32_t face1;             /* Second intersecting face */
    geomlib_vec3 point;         /* Approximate intersection point */
} geomlib_intersection_pair;

/* Self-intersection result */
typedef struct geomlib_self_intersection {
    bool has_intersections;
    geomlib_intersection_pair *pairs;
    size_t pair_count;
} geomlib_self_intersection;

/* Detect self-intersections */
geomlib_self_intersection geomlib_detect_self_intersections(
    const geomlib_trimesh *mesh
);

/* Free self-intersection result */
void geomlib_self_intersection_free(geomlib_self_intersection *result);

/* ============================================================================
 * Adjacency Queries
 * ============================================================================ */

/* Vertex adjacency (one-ring neighborhood) */
typedef struct geomlib_vertex_adjacency {
    uint32_t **neighbors;       /* Neighbor vertices for each vertex */
    uint32_t *neighbor_counts;

    uint32_t **faces;           /* Adjacent faces for each vertex */
    uint32_t *face_counts;
} geomlib_vertex_adjacency;

/* Build vertex adjacency */
geomlib_vertex_adjacency geomlib_build_vertex_adjacency(
    const geomlib_trimesh *mesh
);

/* Free vertex adjacency */
void geomlib_vertex_adjacency_free(geomlib_vertex_adjacency *adjacency);

/* Face adjacency (edge neighbors) */
typedef struct geomlib_face_adjacency {
    uint32_t *neighbors;        /* 3 neighbors per face (UINT32_MAX if boundary) */
} geomlib_face_adjacency;

/* Build face adjacency */
geomlib_face_adjacency geomlib_build_face_adjacency(
    const geomlib_trimesh *mesh
);

/* Free face adjacency */
void geomlib_face_adjacency_free(geomlib_face_adjacency *adjacency);

/* ============================================================================
 * Mesh Orientation
 * ============================================================================ */

/* Check if mesh has consistent orientation */
bool geomlib_is_consistently_oriented(const geomlib_trimesh *mesh);

/* Make mesh consistently oriented (modifies mesh) */
int geomlib_orient_consistently(geomlib_trimesh *mesh);

/* Reverse mesh orientation (flip all normals) */
int geomlib_reverse_orientation(geomlib_trimesh *mesh);

/* ============================================================================
 * Topological Operations
 * ============================================================================ */

/* Remove isolated vertices (not referenced by any face) */
int geomlib_remove_isolated_vertices(geomlib_trimesh *mesh);

/* Remove duplicate faces */
int geomlib_remove_duplicate_faces(geomlib_trimesh *mesh);

/* Split disconnected incident fans. Returns the number of splits or a negative error. */
int geomlib_split_nonmanifold_vertices(geomlib_trimesh *mesh);

/* ============================================================================
 * Valence Analysis
 * ============================================================================ */

/* Vertex valence statistics */
typedef struct geomlib_valence_stats {
    uint32_t *valences;         /* Valence for each vertex */
    uint32_t min_valence;
    uint32_t max_valence;
    float avg_valence;

    /* Histogram */
    uint32_t *histogram;        /* histogram[i] = count of vertices with valence i */
    size_t histogram_size;
} geomlib_valence_stats;

/* Compute vertex valences */
geomlib_valence_stats geomlib_compute_valence_stats(
    const geomlib_trimesh *mesh
);

/* Free valence stats */
void geomlib_valence_stats_free(geomlib_valence_stats *stats);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_TOPOLOGY_H */
