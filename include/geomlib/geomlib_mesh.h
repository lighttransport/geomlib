/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Mesh and Polygon Data Structures
 * Core functions implemented in C for maximum portability
 */

#ifndef GEOMLIB_MESH_H
#define GEOMLIB_MESH_H

#include <stddef.h>
#include <stdint.h>
#include "geomlib_linalg.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Return Codes
 * ============================================================================ */

#define GEOMLIB_MESH_SUCCESS          0
#define GEOMLIB_MESH_ERROR_MEMORY     1
#define GEOMLIB_MESH_ERROR_INVALID    2
#define GEOMLIB_MESH_ERROR_DEGENERATE 3

/* ============================================================================
 * Basic Polygon Types
 * ============================================================================ */

/**
 * Simple 2D polygon (vertex list)
 */
typedef struct geomlib_polygon2d {
    geomlib_vec2 *vertices;
    size_t vertex_count;
    size_t capacity;
} geomlib_polygon2d;

/**
 * Simple 3D polygon (vertex list)
 */
typedef struct geomlib_polygon3d {
    geomlib_vec3 *vertices;
    size_t vertex_count;
    size_t capacity;
} geomlib_polygon3d;

/* ============================================================================
 * Triangle Mesh
 * ============================================================================ */

/**
 * Triangle mesh with indexed vertices
 */
typedef struct geomlib_trimesh {
    geomlib_vec3 *vertices;     /* Vertex positions */
    geomlib_vec3 *normals;      /* Vertex normals (optional, can be NULL) */
    geomlib_vec2 *texcoords;    /* Texture coordinates (optional, can be NULL) */
    uint32_t *indices;          /* Triangle indices (3 per triangle) */
    
    size_t vertex_count;
    size_t triangle_count;
    size_t vertex_capacity;
    size_t triangle_capacity;
} geomlib_trimesh;

/**
 * Double-precision triangle mesh
 */
typedef struct geomlib_trimesh_d {
    geomlib_vec3d *vertices;
    geomlib_vec3d *normals;
    geomlib_vec2d *texcoords;
    uint32_t *indices;
    
    size_t vertex_count;
    size_t triangle_count;
    size_t vertex_capacity;
    size_t triangle_capacity;
} geomlib_trimesh_d;

/* ============================================================================
 * Polygon Mesh (General n-gon support)
 * ============================================================================ */

/**
 * Face structure for polygon mesh
 * Each face can have variable number of vertices
 */
typedef struct geomlib_face {
    uint32_t *vertex_indices;   /* Vertex indices for this face */
    size_t vertex_count;        /* Number of vertices in this face */
} geomlib_face;

/**
 * General polygon mesh supporting n-gons
 */
typedef struct geomlib_polymesh {
    geomlib_vec3 *vertices;
    geomlib_vec3 *normals;      /* Per-vertex normals (optional) */
    geomlib_vec2 *texcoords;    /* Per-vertex texcoords (optional) */
    
    geomlib_face *faces;
    
    size_t vertex_count;
    size_t face_count;
    size_t vertex_capacity;
    size_t face_capacity;
} geomlib_polymesh;

/* ============================================================================
 * Axis-Aligned Bounding Box
 * ============================================================================ */

typedef struct geomlib_aabb {
    geomlib_vec3 min;
    geomlib_vec3 max;
} geomlib_aabb;

typedef struct geomlib_aabb_d {
    geomlib_vec3d min;
    geomlib_vec3d max;
} geomlib_aabb_d;

/* ============================================================================
 * Polygon2D Operations
 * ============================================================================ */

int geomlib_polygon2d_init(geomlib_polygon2d *poly, size_t initial_capacity);
void geomlib_polygon2d_free(geomlib_polygon2d *poly);
int geomlib_polygon2d_add_vertex(geomlib_polygon2d *poly, geomlib_vec2 v);
int geomlib_polygon2d_clear(geomlib_polygon2d *poly);
float geomlib_polygon2d_area(const geomlib_polygon2d *poly);
float geomlib_polygon2d_perimeter(const geomlib_polygon2d *poly);
geomlib_vec2 geomlib_polygon2d_centroid(const geomlib_polygon2d *poly);
int geomlib_polygon2d_is_convex(const geomlib_polygon2d *poly);
int geomlib_polygon2d_is_clockwise(const geomlib_polygon2d *poly);
int geomlib_polygon2d_contains_point(const geomlib_polygon2d *poly, geomlib_vec2 point);

/* ============================================================================
 * Polygon3D Operations
 * ============================================================================ */

int geomlib_polygon3d_init(geomlib_polygon3d *poly, size_t initial_capacity);
void geomlib_polygon3d_free(geomlib_polygon3d *poly);
int geomlib_polygon3d_add_vertex(geomlib_polygon3d *poly, geomlib_vec3 v);
int geomlib_polygon3d_clear(geomlib_polygon3d *poly);
float geomlib_polygon3d_area(const geomlib_polygon3d *poly);
geomlib_vec3 geomlib_polygon3d_normal(const geomlib_polygon3d *poly);
geomlib_vec3 geomlib_polygon3d_centroid(const geomlib_polygon3d *poly);

/* ============================================================================
 * Triangle Mesh Operations
 * ============================================================================ */

int geomlib_trimesh_init(geomlib_trimesh *mesh, size_t vertex_capacity, size_t triangle_capacity);
void geomlib_trimesh_free(geomlib_trimesh *mesh);
int geomlib_trimesh_add_vertex(geomlib_trimesh *mesh, geomlib_vec3 pos, geomlib_vec3 *normal, geomlib_vec2 *texcoord);
int geomlib_trimesh_add_triangle(geomlib_trimesh *mesh, uint32_t i0, uint32_t i1, uint32_t i2);
int geomlib_trimesh_reserve_vertices(geomlib_trimesh *mesh, size_t count);
int geomlib_trimesh_reserve_triangles(geomlib_trimesh *mesh, size_t count);
void geomlib_trimesh_clear(geomlib_trimesh *mesh);
int geomlib_trimesh_compute_normals(geomlib_trimesh *mesh);
geomlib_aabb geomlib_trimesh_compute_aabb(const geomlib_trimesh *mesh);
float geomlib_trimesh_surface_area(const geomlib_trimesh *mesh);
float geomlib_trimesh_volume(const geomlib_trimesh *mesh);
geomlib_vec3 geomlib_trimesh_centroid(const geomlib_trimesh *mesh);

/* Triangle access helpers */
void geomlib_trimesh_get_triangle_vertices(const geomlib_trimesh *mesh, size_t tri_index,
                                           geomlib_vec3 *v0, geomlib_vec3 *v1, geomlib_vec3 *v2);
geomlib_vec3 geomlib_trimesh_triangle_normal(const geomlib_trimesh *mesh, size_t tri_index);
float geomlib_trimesh_triangle_area(const geomlib_trimesh *mesh, size_t tri_index);

/* ============================================================================
 * Polygon Mesh Operations
 * ============================================================================ */

int geomlib_polymesh_init(geomlib_polymesh *mesh, size_t vertex_capacity, size_t face_capacity);
void geomlib_polymesh_free(geomlib_polymesh *mesh);
int geomlib_polymesh_add_vertex(geomlib_polymesh *mesh, geomlib_vec3 pos);
int geomlib_polymesh_add_face(geomlib_polymesh *mesh, const uint32_t *indices, size_t count);
void geomlib_polymesh_clear(geomlib_polymesh *mesh);

/* Conversion */
int geomlib_polymesh_triangulate(const geomlib_polymesh *poly_mesh, geomlib_trimesh *tri_mesh);

/* ============================================================================
 * AABB Operations
 * ============================================================================ */

geomlib_aabb geomlib_aabb_create(geomlib_vec3 min, geomlib_vec3 max);
geomlib_aabb geomlib_aabb_empty(void);
geomlib_aabb geomlib_aabb_expand(geomlib_aabb box, geomlib_vec3 point);
geomlib_aabb geomlib_aabb_merge(geomlib_aabb a, geomlib_aabb b);
int geomlib_aabb_contains_point(geomlib_aabb box, geomlib_vec3 point);
int geomlib_aabb_intersects(geomlib_aabb a, geomlib_aabb b);
geomlib_vec3 geomlib_aabb_center(geomlib_aabb box);
geomlib_vec3 geomlib_aabb_size(geomlib_aabb box);
float geomlib_aabb_surface_area(geomlib_aabb box);
float geomlib_aabb_volume(geomlib_aabb box);
geomlib_vec3 geomlib_aabb_closest_point(const geomlib_aabb *box, geomlib_vec3 point);

/* ============================================================================
 * Primitive Generation
 * ============================================================================ */

int geomlib_create_box(geomlib_trimesh *mesh, float width, float height, float depth);
int geomlib_create_sphere(geomlib_trimesh *mesh, float radius, int segments, int rings);
int geomlib_create_cylinder(geomlib_trimesh *mesh, float radius, float height, int segments);
int geomlib_create_plane(geomlib_trimesh *mesh, float width, float height, int subdivisions);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_MESH_H */
