/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Mesh and Polygon Implementation
 */

#include "geomlib/geomlib_mesh.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Polygon2D Operations
 * ============================================================================ */

int geomlib_polygon2d_init(geomlib_polygon2d *poly, size_t initial_capacity) {
    if (poly == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (initial_capacity == 0) {
        initial_capacity = 8;
    }
    
    poly->vertices = (geomlib_vec2 *)malloc(initial_capacity * sizeof(geomlib_vec2));
    if (poly->vertices == NULL) {
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    
    poly->vertex_count = 0;
    poly->capacity = initial_capacity;
    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_polygon2d_free(geomlib_polygon2d *poly) {
    if (poly != NULL) {
        free(poly->vertices);
        poly->vertices = NULL;
        poly->vertex_count = 0;
        poly->capacity = 0;
    }
}

int geomlib_polygon2d_add_vertex(geomlib_polygon2d *poly, geomlib_vec2 v) {
    if (poly == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (poly->vertex_count >= poly->capacity) {
        size_t new_capacity = poly->capacity * 2;
        geomlib_vec2 *new_vertices = (geomlib_vec2 *)realloc(poly->vertices, 
                                                              new_capacity * sizeof(geomlib_vec2));
        if (new_vertices == NULL) {
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        poly->vertices = new_vertices;
        poly->capacity = new_capacity;
    }
    
    poly->vertices[poly->vertex_count++] = v;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_polygon2d_clear(geomlib_polygon2d *poly) {
    if (poly == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    poly->vertex_count = 0;
    return GEOMLIB_MESH_SUCCESS;
}

float geomlib_polygon2d_area(const geomlib_polygon2d *poly) {
    if (poly == NULL || poly->vertex_count < 3) {
        return 0.0f;
    }
    
    /* Shoelace formula */
    float area = 0.0f;
    size_t n = poly->vertex_count;
    
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        area += poly->vertices[i].x * poly->vertices[j].y;
        area -= poly->vertices[j].x * poly->vertices[i].y;
    }
    
    return fabsf(area) * 0.5f;
}

float geomlib_polygon2d_perimeter(const geomlib_polygon2d *poly) {
    if (poly == NULL || poly->vertex_count < 2) {
        return 0.0f;
    }
    
    float perimeter = 0.0f;
    size_t n = poly->vertex_count;
    
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        geomlib_vec2 edge = geomlib_vec2_sub(poly->vertices[j], poly->vertices[i]);
        perimeter += geomlib_vec2_length(edge);
    }
    
    return perimeter;
}

geomlib_vec2 geomlib_polygon2d_centroid(const geomlib_polygon2d *poly) {
    geomlib_vec2 centroid = {0.0f, 0.0f};
    
    if (poly == NULL || poly->vertex_count == 0) {
        return centroid;
    }
    
    for (size_t i = 0; i < poly->vertex_count; i++) {
        centroid.x += poly->vertices[i].x;
        centroid.y += poly->vertices[i].y;
    }
    
    float inv_count = 1.0f / (float)poly->vertex_count;
    centroid.x *= inv_count;
    centroid.y *= inv_count;
    
    return centroid;
}

int geomlib_polygon2d_is_convex(const geomlib_polygon2d *poly) {
    if (poly == NULL || poly->vertex_count < 3) {
        return 0;
    }
    
    size_t n = poly->vertex_count;
    int sign = 0;
    
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        size_t k = (i + 2) % n;
        
        geomlib_vec2 d1 = geomlib_vec2_sub(poly->vertices[j], poly->vertices[i]);
        geomlib_vec2 d2 = geomlib_vec2_sub(poly->vertices[k], poly->vertices[j]);
        
        float cross = d1.x * d2.y - d1.y * d2.x;
        
        if (fabsf(cross) > 1e-10f) {
            int new_sign = (cross > 0) ? 1 : -1;
            if (sign == 0) {
                sign = new_sign;
            } else if (sign != new_sign) {
                return 0;
            }
        }
    }
    
    return 1;
}

int geomlib_polygon2d_is_clockwise(const geomlib_polygon2d *poly) {
    if (poly == NULL || poly->vertex_count < 3) {
        return 0;
    }
    
    /* Use signed area - negative means clockwise */
    float signed_area = 0.0f;
    size_t n = poly->vertex_count;
    
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        signed_area += poly->vertices[i].x * poly->vertices[j].y;
        signed_area -= poly->vertices[j].x * poly->vertices[i].y;
    }
    
    return signed_area < 0.0f;
}

int geomlib_polygon2d_contains_point(const geomlib_polygon2d *poly, geomlib_vec2 point) {
    if (poly == NULL || poly->vertex_count < 3) {
        return 0;
    }
    
    /* Ray casting algorithm */
    int inside = 0;
    size_t n = poly->vertex_count;
    
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        float yi = poly->vertices[i].y;
        float yj = poly->vertices[j].y;
        
        if ((yi > point.y) != (yj > point.y)) {
            float xi = poly->vertices[i].x;
            float xj = poly->vertices[j].x;
            float x_intersect = (xj - xi) * (point.y - yi) / (yj - yi) + xi;
            
            if (point.x < x_intersect) {
                inside = !inside;
            }
        }
    }
    
    return inside;
}

/* ============================================================================
 * Polygon3D Operations
 * ============================================================================ */

int geomlib_polygon3d_init(geomlib_polygon3d *poly, size_t initial_capacity) {
    if (poly == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (initial_capacity == 0) {
        initial_capacity = 8;
    }
    
    poly->vertices = (geomlib_vec3 *)malloc(initial_capacity * sizeof(geomlib_vec3));
    if (poly->vertices == NULL) {
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    
    poly->vertex_count = 0;
    poly->capacity = initial_capacity;
    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_polygon3d_free(geomlib_polygon3d *poly) {
    if (poly != NULL) {
        free(poly->vertices);
        poly->vertices = NULL;
        poly->vertex_count = 0;
        poly->capacity = 0;
    }
}

int geomlib_polygon3d_add_vertex(geomlib_polygon3d *poly, geomlib_vec3 v) {
    if (poly == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (poly->vertex_count >= poly->capacity) {
        size_t new_capacity = poly->capacity * 2;
        geomlib_vec3 *new_vertices = (geomlib_vec3 *)realloc(poly->vertices, 
                                                              new_capacity * sizeof(geomlib_vec3));
        if (new_vertices == NULL) {
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        poly->vertices = new_vertices;
        poly->capacity = new_capacity;
    }
    
    poly->vertices[poly->vertex_count++] = v;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_polygon3d_clear(geomlib_polygon3d *poly) {
    if (poly == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    poly->vertex_count = 0;
    return GEOMLIB_MESH_SUCCESS;
}

float geomlib_polygon3d_area(const geomlib_polygon3d *poly) {
    if (poly == NULL || poly->vertex_count < 3) {
        return 0.0f;
    }
    
    /* Compute area using cross products */
    geomlib_vec3 sum = {0.0f, 0.0f, 0.0f};
    size_t n = poly->vertex_count;
    
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        geomlib_vec3 cross = geomlib_vec3_cross(poly->vertices[i], poly->vertices[j]);
        sum = geomlib_vec3_add(sum, cross);
    }
    
    return geomlib_vec3_length(sum) * 0.5f;
}

geomlib_vec3 geomlib_polygon3d_normal(const geomlib_polygon3d *poly) {
    geomlib_vec3 normal = {0.0f, 0.0f, 0.0f};
    
    if (poly == NULL || poly->vertex_count < 3) {
        return normal;
    }
    
    /* Newell's method for computing polygon normal */
    size_t n = poly->vertex_count;
    
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        normal.x += (poly->vertices[i].y - poly->vertices[j].y) * 
                    (poly->vertices[i].z + poly->vertices[j].z);
        normal.y += (poly->vertices[i].z - poly->vertices[j].z) * 
                    (poly->vertices[i].x + poly->vertices[j].x);
        normal.z += (poly->vertices[i].x - poly->vertices[j].x) * 
                    (poly->vertices[i].y + poly->vertices[j].y);
    }
    
    return geomlib_vec3_normalize(normal);
}

geomlib_vec3 geomlib_polygon3d_centroid(const geomlib_polygon3d *poly) {
    geomlib_vec3 centroid = {0.0f, 0.0f, 0.0f};
    
    if (poly == NULL || poly->vertex_count == 0) {
        return centroid;
    }
    
    for (size_t i = 0; i < poly->vertex_count; i++) {
        centroid = geomlib_vec3_add(centroid, poly->vertices[i]);
    }
    
    return geomlib_vec3_div(centroid, (float)poly->vertex_count);
}

/* ============================================================================
 * Triangle Mesh Operations
 * ============================================================================ */

int geomlib_trimesh_init(geomlib_trimesh *mesh, size_t vertex_capacity, size_t triangle_capacity) {
    if (mesh == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (vertex_capacity == 0) vertex_capacity = 64;
    if (triangle_capacity == 0) triangle_capacity = 64;
    
    mesh->vertices = (geomlib_vec3 *)malloc(vertex_capacity * sizeof(geomlib_vec3));
    mesh->normals = NULL;  /* Allocated on demand */
    mesh->texcoords = NULL;
    mesh->indices = (uint32_t *)malloc(triangle_capacity * 3 * sizeof(uint32_t));
    
    if (mesh->vertices == NULL || mesh->indices == NULL) {
        free(mesh->vertices);
        free(mesh->indices);
        mesh->vertices = NULL;
        mesh->indices = NULL;
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    
    mesh->vertex_count = 0;
    mesh->triangle_count = 0;
    mesh->vertex_capacity = vertex_capacity;
    mesh->triangle_capacity = triangle_capacity;
    
    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_trimesh_free(geomlib_trimesh *mesh) {
    if (mesh != NULL) {
        free(mesh->vertices);
        free(mesh->normals);
        free(mesh->texcoords);
        free(mesh->indices);
        mesh->vertices = NULL;
        mesh->normals = NULL;
        mesh->texcoords = NULL;
        mesh->indices = NULL;
        mesh->vertex_count = 0;
        mesh->triangle_count = 0;
        mesh->vertex_capacity = 0;
        mesh->triangle_capacity = 0;
    }
}

int geomlib_trimesh_reserve_vertices(geomlib_trimesh *mesh, size_t count) {
    if (mesh == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (count > mesh->vertex_capacity) {
        geomlib_vec3 *new_vertices = (geomlib_vec3 *)realloc(mesh->vertices, 
                                                              count * sizeof(geomlib_vec3));
        if (new_vertices == NULL) {
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        mesh->vertices = new_vertices;
        
        if (mesh->normals != NULL) {
            geomlib_vec3 *new_normals = (geomlib_vec3 *)realloc(mesh->normals,
                                                                 count * sizeof(geomlib_vec3));
            if (new_normals == NULL) {
                return GEOMLIB_MESH_ERROR_MEMORY;
            }
            mesh->normals = new_normals;
        }
        
        if (mesh->texcoords != NULL) {
            geomlib_vec2 *new_texcoords = (geomlib_vec2 *)realloc(mesh->texcoords,
                                                                   count * sizeof(geomlib_vec2));
            if (new_texcoords == NULL) {
                return GEOMLIB_MESH_ERROR_MEMORY;
            }
            mesh->texcoords = new_texcoords;
        }
        
        mesh->vertex_capacity = count;
    }
    
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_trimesh_reserve_triangles(geomlib_trimesh *mesh, size_t count) {
    if (mesh == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (count > mesh->triangle_capacity) {
        uint32_t *new_indices = (uint32_t *)realloc(mesh->indices, 
                                                     count * 3 * sizeof(uint32_t));
        if (new_indices == NULL) {
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        mesh->indices = new_indices;
        mesh->triangle_capacity = count;
    }
    
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_trimesh_add_vertex(geomlib_trimesh *mesh, geomlib_vec3 pos, 
                               geomlib_vec3 *normal, geomlib_vec2 *texcoord) {
    if (mesh == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (mesh->vertex_count >= mesh->vertex_capacity) {
        int result = geomlib_trimesh_reserve_vertices(mesh, mesh->vertex_capacity * 2);
        if (result != GEOMLIB_MESH_SUCCESS) {
            return result;
        }
    }
    
    mesh->vertices[mesh->vertex_count] = pos;
    
    if (normal != NULL) {
        if (mesh->normals == NULL) {
            mesh->normals = (geomlib_vec3 *)calloc(mesh->vertex_capacity, sizeof(geomlib_vec3));
            if (mesh->normals == NULL) {
                return GEOMLIB_MESH_ERROR_MEMORY;
            }
        }
        mesh->normals[mesh->vertex_count] = *normal;
    }
    
    if (texcoord != NULL) {
        if (mesh->texcoords == NULL) {
            mesh->texcoords = (geomlib_vec2 *)calloc(mesh->vertex_capacity, sizeof(geomlib_vec2));
            if (mesh->texcoords == NULL) {
                return GEOMLIB_MESH_ERROR_MEMORY;
            }
        }
        mesh->texcoords[mesh->vertex_count] = *texcoord;
    }
    
    mesh->vertex_count++;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_trimesh_add_triangle(geomlib_trimesh *mesh, uint32_t i0, uint32_t i1, uint32_t i2) {
    if (mesh == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (mesh->triangle_count >= mesh->triangle_capacity) {
        int result = geomlib_trimesh_reserve_triangles(mesh, mesh->triangle_capacity * 2);
        if (result != GEOMLIB_MESH_SUCCESS) {
            return result;
        }
    }
    
    size_t base = mesh->triangle_count * 3;
    mesh->indices[base + 0] = i0;
    mesh->indices[base + 1] = i1;
    mesh->indices[base + 2] = i2;
    mesh->triangle_count++;
    
    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_trimesh_clear(geomlib_trimesh *mesh) {
    if (mesh != NULL) {
        mesh->vertex_count = 0;
        mesh->triangle_count = 0;
    }
}

void geomlib_trimesh_get_triangle_vertices(const geomlib_trimesh *mesh, size_t tri_index,
                                           geomlib_vec3 *v0, geomlib_vec3 *v1, geomlib_vec3 *v2) {
    if (mesh == NULL || tri_index >= mesh->triangle_count) {
        return;
    }
    
    size_t base = tri_index * 3;
    if (v0) *v0 = mesh->vertices[mesh->indices[base + 0]];
    if (v1) *v1 = mesh->vertices[mesh->indices[base + 1]];
    if (v2) *v2 = mesh->vertices[mesh->indices[base + 2]];
}

geomlib_vec3 geomlib_trimesh_triangle_normal(const geomlib_trimesh *mesh, size_t tri_index) {
    geomlib_vec3 normal = {0.0f, 0.0f, 0.0f};
    
    if (mesh == NULL || tri_index >= mesh->triangle_count) {
        return normal;
    }
    
    geomlib_vec3 v0, v1, v2;
    geomlib_trimesh_get_triangle_vertices(mesh, tri_index, &v0, &v1, &v2);
    
    geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
    
    return geomlib_vec3_normalize(geomlib_vec3_cross(e1, e2));
}

float geomlib_trimesh_triangle_area(const geomlib_trimesh *mesh, size_t tri_index) {
    if (mesh == NULL || tri_index >= mesh->triangle_count) {
        return 0.0f;
    }
    
    geomlib_vec3 v0, v1, v2;
    geomlib_trimesh_get_triangle_vertices(mesh, tri_index, &v0, &v1, &v2);
    
    geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
    geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
    geomlib_vec3 cross = geomlib_vec3_cross(e1, e2);
    
    return geomlib_vec3_length(cross) * 0.5f;
}

int geomlib_trimesh_compute_normals(geomlib_trimesh *mesh) {
    if (mesh == NULL || mesh->vertex_count == 0 || mesh->triangle_count == 0) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    /* Allocate normals if needed */
    if (mesh->normals == NULL) {
        mesh->normals = (geomlib_vec3 *)calloc(mesh->vertex_capacity, sizeof(geomlib_vec3));
        if (mesh->normals == NULL) {
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
    }
    
    /* Zero out normals */
    memset(mesh->normals, 0, mesh->vertex_count * sizeof(geomlib_vec3));
    
    /* Accumulate face normals weighted by area */
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        size_t base = i * 3;
        uint32_t i0 = mesh->indices[base + 0];
        uint32_t i1 = mesh->indices[base + 1];
        uint32_t i2 = mesh->indices[base + 2];
        
        geomlib_vec3 v0 = mesh->vertices[i0];
        geomlib_vec3 v1 = mesh->vertices[i1];
        geomlib_vec3 v2 = mesh->vertices[i2];
        
        geomlib_vec3 e1 = geomlib_vec3_sub(v1, v0);
        geomlib_vec3 e2 = geomlib_vec3_sub(v2, v0);
        geomlib_vec3 face_normal = geomlib_vec3_cross(e1, e2);
        
        /* Add unnormalized face normal (weighted by area) */
        mesh->normals[i0] = geomlib_vec3_add(mesh->normals[i0], face_normal);
        mesh->normals[i1] = geomlib_vec3_add(mesh->normals[i1], face_normal);
        mesh->normals[i2] = geomlib_vec3_add(mesh->normals[i2], face_normal);
    }
    
    /* Normalize vertex normals */
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        mesh->normals[i] = geomlib_vec3_normalize(mesh->normals[i]);
    }
    
    return GEOMLIB_MESH_SUCCESS;
}

geomlib_aabb geomlib_trimesh_compute_aabb(const geomlib_trimesh *mesh) {
    geomlib_aabb box = geomlib_aabb_empty();
    
    if (mesh == NULL || mesh->vertex_count == 0) {
        return box;
    }
    
    box.min = mesh->vertices[0];
    box.max = mesh->vertices[0];
    
    for (size_t i = 1; i < mesh->vertex_count; i++) {
        box = geomlib_aabb_expand(box, mesh->vertices[i]);
    }
    
    return box;
}

float geomlib_trimesh_surface_area(const geomlib_trimesh *mesh) {
    if (mesh == NULL) {
        return 0.0f;
    }
    
    float total_area = 0.0f;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        total_area += geomlib_trimesh_triangle_area(mesh, i);
    }
    
    return total_area;
}

float geomlib_trimesh_volume(const geomlib_trimesh *mesh) {
    if (mesh == NULL || mesh->triangle_count == 0) {
        return 0.0f;
    }
    
    /* Compute volume using signed tetrahedron volumes */
    float volume = 0.0f;
    
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 v0, v1, v2;
        geomlib_trimesh_get_triangle_vertices(mesh, i, &v0, &v1, &v2);
        
        /* Volume of tetrahedron with origin */
        float v = geomlib_vec3_dot(v0, geomlib_vec3_cross(v1, v2));
        volume += v;
    }
    
    return fabsf(volume) / 6.0f;
}

geomlib_vec3 geomlib_trimesh_centroid(const geomlib_trimesh *mesh) {
    geomlib_vec3 centroid = {0.0f, 0.0f, 0.0f};
    
    if (mesh == NULL || mesh->vertex_count == 0) {
        return centroid;
    }
    
    for (size_t i = 0; i < mesh->vertex_count; i++) {
        centroid = geomlib_vec3_add(centroid, mesh->vertices[i]);
    }
    
    return geomlib_vec3_div(centroid, (float)mesh->vertex_count);
}

/* ============================================================================
 * Polygon Mesh Operations
 * ============================================================================ */

int geomlib_polymesh_init(geomlib_polymesh *mesh, size_t vertex_capacity, size_t face_capacity) {
    if (mesh == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (vertex_capacity == 0) vertex_capacity = 64;
    if (face_capacity == 0) face_capacity = 32;
    
    mesh->vertices = (geomlib_vec3 *)malloc(vertex_capacity * sizeof(geomlib_vec3));
    mesh->faces = (geomlib_face *)calloc(face_capacity, sizeof(geomlib_face));
    mesh->normals = NULL;
    mesh->texcoords = NULL;
    
    if (mesh->vertices == NULL || mesh->faces == NULL) {
        free(mesh->vertices);
        free(mesh->faces);
        mesh->vertices = NULL;
        mesh->faces = NULL;
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    
    mesh->vertex_count = 0;
    mesh->face_count = 0;
    mesh->vertex_capacity = vertex_capacity;
    mesh->face_capacity = face_capacity;
    
    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_polymesh_free(geomlib_polymesh *mesh) {
    if (mesh != NULL) {
        free(mesh->vertices);
        free(mesh->normals);
        free(mesh->texcoords);
        
        if (mesh->faces != NULL) {
            for (size_t i = 0; i < mesh->face_count; i++) {
                free(mesh->faces[i].vertex_indices);
            }
            free(mesh->faces);
        }
        
        mesh->vertices = NULL;
        mesh->normals = NULL;
        mesh->texcoords = NULL;
        mesh->faces = NULL;
        mesh->vertex_count = 0;
        mesh->face_count = 0;
    }
}

int geomlib_polymesh_add_vertex(geomlib_polymesh *mesh, geomlib_vec3 pos) {
    if (mesh == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (mesh->vertex_count >= mesh->vertex_capacity) {
        size_t new_capacity = mesh->vertex_capacity * 2;
        geomlib_vec3 *new_vertices = (geomlib_vec3 *)realloc(mesh->vertices,
                                                              new_capacity * sizeof(geomlib_vec3));
        if (new_vertices == NULL) {
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        mesh->vertices = new_vertices;
        mesh->vertex_capacity = new_capacity;
    }
    
    mesh->vertices[mesh->vertex_count++] = pos;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_polymesh_add_face(geomlib_polymesh *mesh, const uint32_t *indices, size_t count) {
    if (mesh == NULL || indices == NULL || count < 3) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    if (mesh->face_count >= mesh->face_capacity) {
        size_t new_capacity = mesh->face_capacity * 2;
        geomlib_face *new_faces = (geomlib_face *)realloc(mesh->faces,
                                                           new_capacity * sizeof(geomlib_face));
        if (new_faces == NULL) {
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        mesh->faces = new_faces;
        mesh->face_capacity = new_capacity;
        
        /* Initialize new faces */
        for (size_t i = mesh->face_count; i < new_capacity; i++) {
            mesh->faces[i].vertex_indices = NULL;
            mesh->faces[i].vertex_count = 0;
        }
    }
    
    geomlib_face *face = &mesh->faces[mesh->face_count];
    face->vertex_indices = (uint32_t *)malloc(count * sizeof(uint32_t));
    if (face->vertex_indices == NULL) {
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    
    memcpy(face->vertex_indices, indices, count * sizeof(uint32_t));
    face->vertex_count = count;
    mesh->face_count++;
    
    return GEOMLIB_MESH_SUCCESS;
}

void geomlib_polymesh_clear(geomlib_polymesh *mesh) {
    if (mesh != NULL) {
        for (size_t i = 0; i < mesh->face_count; i++) {
            free(mesh->faces[i].vertex_indices);
            mesh->faces[i].vertex_indices = NULL;
            mesh->faces[i].vertex_count = 0;
        }
        mesh->vertex_count = 0;
        mesh->face_count = 0;
    }
}

int geomlib_polymesh_triangulate(const geomlib_polymesh *poly_mesh, geomlib_trimesh *tri_mesh) {
    if (poly_mesh == NULL || tri_mesh == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    /* Count total triangles needed */
    size_t total_triangles = 0;
    for (size_t i = 0; i < poly_mesh->face_count; i++) {
        if (poly_mesh->faces[i].vertex_count >= 3) {
            total_triangles += poly_mesh->faces[i].vertex_count - 2;
        }
    }
    
    int result = geomlib_trimesh_init(tri_mesh, poly_mesh->vertex_count, total_triangles);
    if (result != GEOMLIB_MESH_SUCCESS) {
        return result;
    }
    
    /* Copy vertices */
    for (size_t i = 0; i < poly_mesh->vertex_count; i++) {
        result = geomlib_trimesh_add_vertex(tri_mesh, poly_mesh->vertices[i], NULL, NULL);
        if (result != GEOMLIB_MESH_SUCCESS) {
            geomlib_trimesh_free(tri_mesh);
            return result;
        }
    }
    
    /* Triangulate faces using fan triangulation */
    for (size_t i = 0; i < poly_mesh->face_count; i++) {
        const geomlib_face *face = &poly_mesh->faces[i];
        if (face->vertex_count < 3) continue;
        
        uint32_t i0 = face->vertex_indices[0];
        for (size_t j = 1; j < face->vertex_count - 1; j++) {
            uint32_t i1 = face->vertex_indices[j];
            uint32_t i2 = face->vertex_indices[j + 1];
            
            result = geomlib_trimesh_add_triangle(tri_mesh, i0, i1, i2);
            if (result != GEOMLIB_MESH_SUCCESS) {
                geomlib_trimesh_free(tri_mesh);
                return result;
            }
        }
    }
    
    return GEOMLIB_MESH_SUCCESS;
}

/* ============================================================================
 * AABB Operations
 * ============================================================================ */

geomlib_aabb geomlib_aabb_create(geomlib_vec3 min, geomlib_vec3 max) {
    geomlib_aabb box;
    box.min = min;
    box.max = max;
    return box;
}

geomlib_aabb geomlib_aabb_empty(void) {
    geomlib_aabb box;
    box.min = geomlib_vec3_create(FLT_MAX, FLT_MAX, FLT_MAX);
    box.max = geomlib_vec3_create(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    return box;
}

geomlib_aabb geomlib_aabb_expand(geomlib_aabb box, geomlib_vec3 point) {
    if (point.x < box.min.x) box.min.x = point.x;
    if (point.y < box.min.y) box.min.y = point.y;
    if (point.z < box.min.z) box.min.z = point.z;
    if (point.x > box.max.x) box.max.x = point.x;
    if (point.y > box.max.y) box.max.y = point.y;
    if (point.z > box.max.z) box.max.z = point.z;
    return box;
}

geomlib_aabb geomlib_aabb_merge(geomlib_aabb a, geomlib_aabb b) {
    geomlib_aabb box;
    box.min.x = (a.min.x < b.min.x) ? a.min.x : b.min.x;
    box.min.y = (a.min.y < b.min.y) ? a.min.y : b.min.y;
    box.min.z = (a.min.z < b.min.z) ? a.min.z : b.min.z;
    box.max.x = (a.max.x > b.max.x) ? a.max.x : b.max.x;
    box.max.y = (a.max.y > b.max.y) ? a.max.y : b.max.y;
    box.max.z = (a.max.z > b.max.z) ? a.max.z : b.max.z;
    return box;
}

int geomlib_aabb_contains_point(geomlib_aabb box, geomlib_vec3 point) {
    return point.x >= box.min.x && point.x <= box.max.x &&
           point.y >= box.min.y && point.y <= box.max.y &&
           point.z >= box.min.z && point.z <= box.max.z;
}

int geomlib_aabb_intersects(geomlib_aabb a, geomlib_aabb b) {
    return a.min.x <= b.max.x && a.max.x >= b.min.x &&
           a.min.y <= b.max.y && a.max.y >= b.min.y &&
           a.min.z <= b.max.z && a.max.z >= b.min.z;
}

geomlib_vec3 geomlib_aabb_center(geomlib_aabb box) {
    return geomlib_vec3_mul(geomlib_vec3_add(box.min, box.max), 0.5f);
}

geomlib_vec3 geomlib_aabb_size(geomlib_aabb box) {
    return geomlib_vec3_sub(box.max, box.min);
}

float geomlib_aabb_surface_area(geomlib_aabb box) {
    geomlib_vec3 size = geomlib_aabb_size(box);
    return 2.0f * (size.x * size.y + size.y * size.z + size.z * size.x);
}

float geomlib_aabb_volume(geomlib_aabb box) {
    geomlib_vec3 size = geomlib_aabb_size(box);
    return size.x * size.y * size.z;
}

/* ============================================================================
 * Primitive Generation
 * ============================================================================ */

int geomlib_create_box(geomlib_trimesh *mesh, float width, float height, float depth) {
    if (mesh == NULL) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    float hw = width * 0.5f;
    float hh = height * 0.5f;
    float hd = depth * 0.5f;
    
    int result = geomlib_trimesh_init(mesh, 24, 12);
    if (result != GEOMLIB_MESH_SUCCESS) {
        return result;
    }
    
    /* Define vertices for each face (with normals) */
    /* Front face */
    geomlib_vec3 n_front = {0, 0, 1};
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw, -hh, hd), &n_front, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create( hw, -hh, hd), &n_front, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create( hw,  hh, hd), &n_front, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw,  hh, hd), &n_front, NULL);
    
    /* Back face */
    geomlib_vec3 n_back = {0, 0, -1};
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create( hw, -hh, -hd), &n_back, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw, -hh, -hd), &n_back, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw,  hh, -hd), &n_back, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create( hw,  hh, -hd), &n_back, NULL);
    
    /* Top face */
    geomlib_vec3 n_top = {0, 1, 0};
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw, hh,  hd), &n_top, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create( hw, hh,  hd), &n_top, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create( hw, hh, -hd), &n_top, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw, hh, -hd), &n_top, NULL);
    
    /* Bottom face */
    geomlib_vec3 n_bottom = {0, -1, 0};
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw, -hh, -hd), &n_bottom, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create( hw, -hh, -hd), &n_bottom, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create( hw, -hh,  hd), &n_bottom, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw, -hh,  hd), &n_bottom, NULL);
    
    /* Right face */
    geomlib_vec3 n_right = {1, 0, 0};
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(hw, -hh,  hd), &n_right, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(hw, -hh, -hd), &n_right, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(hw,  hh, -hd), &n_right, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(hw,  hh,  hd), &n_right, NULL);
    
    /* Left face */
    geomlib_vec3 n_left = {-1, 0, 0};
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw, -hh, -hd), &n_left, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw, -hh,  hd), &n_left, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw,  hh,  hd), &n_left, NULL);
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(-hw,  hh, -hd), &n_left, NULL);
    
    /* Add triangles for each face */
    for (uint32_t i = 0; i < 6; i++) {
        uint32_t base = i * 4;
        geomlib_trimesh_add_triangle(mesh, base + 0, base + 1, base + 2);
        geomlib_trimesh_add_triangle(mesh, base + 0, base + 2, base + 3);
    }
    
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_create_sphere(geomlib_trimesh *mesh, float radius, int segments, int rings) {
    if (mesh == NULL || segments < 3 || rings < 2) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    size_t vertex_count = (size_t)((segments + 1) * (rings + 1));
    size_t triangle_count = (size_t)(segments * rings * 2);
    
    int result = geomlib_trimesh_init(mesh, vertex_count, triangle_count);
    if (result != GEOMLIB_MESH_SUCCESS) {
        return result;
    }
    
    /* Generate vertices */
    for (int r = 0; r <= rings; r++) {
        float phi = (float)M_PI * (float)r / (float)rings;
        float sin_phi = sinf(phi);
        float cos_phi = cosf(phi);
        
        for (int s = 0; s <= segments; s++) {
            float theta = 2.0f * (float)M_PI * (float)s / (float)segments;
            float sin_theta = sinf(theta);
            float cos_theta = cosf(theta);
            
            geomlib_vec3 pos;
            pos.x = radius * sin_phi * cos_theta;
            pos.y = radius * cos_phi;
            pos.z = radius * sin_phi * sin_theta;
            
            geomlib_vec3 normal = geomlib_vec3_normalize(pos);
            
            geomlib_trimesh_add_vertex(mesh, pos, &normal, NULL);
        }
    }
    
    /* Generate triangles */
    for (int r = 0; r < rings; r++) {
        for (int s = 0; s < segments; s++) {
            uint32_t i0 = (uint32_t)(r * (segments + 1) + s);
            uint32_t i1 = i0 + 1;
            uint32_t i2 = i0 + (uint32_t)(segments + 1);
            uint32_t i3 = i2 + 1;
            
            geomlib_trimesh_add_triangle(mesh, i0, i2, i1);
            geomlib_trimesh_add_triangle(mesh, i1, i2, i3);
        }
    }
    
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_create_cylinder(geomlib_trimesh *mesh, float radius, float height, int segments) {
    if (mesh == NULL || segments < 3) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    float hh = height * 0.5f;
    
    /* Vertices: 2 centers + 2 rings of segments */
    size_t vertex_count = (size_t)(2 + segments * 4);
    /* Triangles: top cap + bottom cap + side */
    size_t triangle_count = (size_t)(segments * 4);
    
    int result = geomlib_trimesh_init(mesh, vertex_count, triangle_count);
    if (result != GEOMLIB_MESH_SUCCESS) {
        return result;
    }
    
    /* Top center vertex */
    geomlib_vec3 n_up = {0, 1, 0};
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(0, hh, 0), &n_up, NULL);
    
    /* Bottom center vertex */
    geomlib_vec3 n_down = {0, -1, 0};
    geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(0, -hh, 0), &n_down, NULL);
    
    /* Top ring vertices */
    for (int i = 0; i < segments; i++) {
        float theta = 2.0f * (float)M_PI * (float)i / (float)segments;
        float x = radius * cosf(theta);
        float z = radius * sinf(theta);
        geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(x, hh, z), &n_up, NULL);
    }
    
    /* Bottom ring vertices */
    for (int i = 0; i < segments; i++) {
        float theta = 2.0f * (float)M_PI * (float)i / (float)segments;
        float x = radius * cosf(theta);
        float z = radius * sinf(theta);
        geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(x, -hh, z), &n_down, NULL);
    }
    
    /* Side vertices (with side normals) */
    for (int i = 0; i < segments; i++) {
        float theta = 2.0f * (float)M_PI * (float)i / (float)segments;
        float x = radius * cosf(theta);
        float z = radius * sinf(theta);
        geomlib_vec3 n_side = geomlib_vec3_normalize(geomlib_vec3_create(x, 0, z));
        geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(x, hh, z), &n_side, NULL);
    }
    
    for (int i = 0; i < segments; i++) {
        float theta = 2.0f * (float)M_PI * (float)i / (float)segments;
        float x = radius * cosf(theta);
        float z = radius * sinf(theta);
        geomlib_vec3 n_side = geomlib_vec3_normalize(geomlib_vec3_create(x, 0, z));
        geomlib_trimesh_add_vertex(mesh, geomlib_vec3_create(x, -hh, z), &n_side, NULL);
    }
    
    /* Top cap triangles */
    for (int i = 0; i < segments; i++) {
        uint32_t next = (uint32_t)((i + 1) % segments);
        geomlib_trimesh_add_triangle(mesh, 0, 2 + (uint32_t)i, 2 + next);
    }
    
    /* Bottom cap triangles */
    for (int i = 0; i < segments; i++) {
        uint32_t next = (uint32_t)((i + 1) % segments);
        geomlib_trimesh_add_triangle(mesh, 1, 2 + (uint32_t)segments + next, 2 + (uint32_t)segments + (uint32_t)i);
    }
    
    /* Side triangles */
    uint32_t side_top_start = 2 + (uint32_t)segments * 2;
    uint32_t side_bottom_start = side_top_start + (uint32_t)segments;
    
    for (int i = 0; i < segments; i++) {
        uint32_t next = (uint32_t)((i + 1) % segments);
        uint32_t t0 = side_top_start + (uint32_t)i;
        uint32_t t1 = side_top_start + next;
        uint32_t b0 = side_bottom_start + (uint32_t)i;
        uint32_t b1 = side_bottom_start + next;
        
        geomlib_trimesh_add_triangle(mesh, t0, b0, t1);
        geomlib_trimesh_add_triangle(mesh, t1, b0, b1);
    }
    
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_create_plane(geomlib_trimesh *mesh, float width, float height, int subdivisions) {
    if (mesh == NULL || subdivisions < 1) {
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    
    int n = subdivisions + 1;
    size_t vertex_count = (size_t)(n * n);
    size_t triangle_count = (size_t)(subdivisions * subdivisions * 2);
    
    int result = geomlib_trimesh_init(mesh, vertex_count, triangle_count);
    if (result != GEOMLIB_MESH_SUCCESS) {
        return result;
    }
    
    float hw = width * 0.5f;
    float hh = height * 0.5f;
    geomlib_vec3 normal = {0, 1, 0};
    
    /* Generate vertices */
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            float u = (float)i / (float)subdivisions;
            float v = (float)j / (float)subdivisions;
            
            geomlib_vec3 pos;
            pos.x = -hw + u * width;
            pos.y = 0;
            pos.z = -hh + v * height;
            
            geomlib_vec2 texcoord = {u, v};
            geomlib_trimesh_add_vertex(mesh, pos, &normal, &texcoord);
        }
    }
    
    /* Generate triangles */
    for (int j = 0; j < subdivisions; j++) {
        for (int i = 0; i < subdivisions; i++) {
            uint32_t i0 = (uint32_t)(j * n + i);
            uint32_t i1 = i0 + 1;
            uint32_t i2 = i0 + (uint32_t)n;
            uint32_t i3 = i2 + 1;
            
            geomlib_trimesh_add_triangle(mesh, i0, i2, i1);
            geomlib_trimesh_add_triangle(mesh, i1, i2, i3);
        }
    }
    
    return GEOMLIB_MESH_SUCCESS;
}
