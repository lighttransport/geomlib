/*
 * geomlib - Mesh Tests
 */

#include <cstdio>
#include <cmath>
#include "geomlib/geomlib_mesh.h"

static int test_count = 0;
static int pass_count = 0;

#define EPSILON 1e-5f

static bool float_eq(float a, float b, float eps = EPSILON) {
    return std::fabs(a - b) < eps;
}

static void test(bool condition, const char* name) {
    test_count++;
    if (condition) {
        pass_count++;
        std::printf("[PASS] %s\n", name);
    } else {
        std::printf("[FAIL] %s\n", name);
    }
}

static void test_polygon2d() {
    geomlib_polygon2d poly;
    int result = geomlib_polygon2d_init(&poly, 4);
    test(result == GEOMLIB_MESH_SUCCESS, "polygon2d_init");
    
    // Create a unit square
    geomlib_polygon2d_add_vertex(&poly, geomlib_vec2_create(0.0f, 0.0f));
    geomlib_polygon2d_add_vertex(&poly, geomlib_vec2_create(1.0f, 0.0f));
    geomlib_polygon2d_add_vertex(&poly, geomlib_vec2_create(1.0f, 1.0f));
    geomlib_polygon2d_add_vertex(&poly, geomlib_vec2_create(0.0f, 1.0f));
    
    test(poly.vertex_count == 4, "polygon2d_vertex_count");
    
    float area = geomlib_polygon2d_area(&poly);
    test(float_eq(area, 1.0f), "polygon2d_area_unit_square");
    
    float perimeter = geomlib_polygon2d_perimeter(&poly);
    test(float_eq(perimeter, 4.0f), "polygon2d_perimeter_unit_square");
    
    geomlib_vec2 centroid = geomlib_polygon2d_centroid(&poly);
    test(float_eq(centroid.x, 0.5f) && float_eq(centroid.y, 0.5f), "polygon2d_centroid");
    
    test(geomlib_polygon2d_is_convex(&poly), "polygon2d_is_convex_square");
    
    // Test point containment
    test(geomlib_polygon2d_contains_point(&poly, geomlib_vec2_create(0.5f, 0.5f)), 
         "polygon2d_contains_inside");
    test(!geomlib_polygon2d_contains_point(&poly, geomlib_vec2_create(2.0f, 2.0f)), 
         "polygon2d_contains_outside");
    
    geomlib_polygon2d_free(&poly);
}

static void test_polygon3d() {
    geomlib_polygon3d poly;
    int result = geomlib_polygon3d_init(&poly, 4);
    test(result == GEOMLIB_MESH_SUCCESS, "polygon3d_init");
    
    // Create a unit square in XZ plane
    geomlib_polygon3d_add_vertex(&poly, geomlib_vec3_create(0.0f, 0.0f, 0.0f));
    geomlib_polygon3d_add_vertex(&poly, geomlib_vec3_create(1.0f, 0.0f, 0.0f));
    geomlib_polygon3d_add_vertex(&poly, geomlib_vec3_create(1.0f, 0.0f, 1.0f));
    geomlib_polygon3d_add_vertex(&poly, geomlib_vec3_create(0.0f, 0.0f, 1.0f));
    
    float area = geomlib_polygon3d_area(&poly);
    test(float_eq(area, 1.0f), "polygon3d_area_unit_square");
    
    geomlib_vec3 normal = geomlib_polygon3d_normal(&poly);
    // Normal should point in Y direction
    test(float_eq(std::fabs(normal.y), 1.0f, 0.01f), "polygon3d_normal");
    
    geomlib_polygon3d_free(&poly);
}

static void test_trimesh_init() {
    geomlib_trimesh mesh;
    int result = geomlib_trimesh_init(&mesh, 64, 64);
    test(result == GEOMLIB_MESH_SUCCESS, "trimesh_init");
    test(mesh.vertex_count == 0, "trimesh_init_vertex_count");
    test(mesh.triangle_count == 0, "trimesh_init_triangle_count");
    geomlib_trimesh_free(&mesh);
}

static void test_trimesh_triangle() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 64, 64);
    
    // Add a simple triangle
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 1.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    
    test(mesh.vertex_count == 3, "trimesh_add_vertices");
    test(mesh.triangle_count == 1, "trimesh_add_triangle");
    
    float area = geomlib_trimesh_triangle_area(&mesh, 0);
    test(float_eq(area, 0.5f), "trimesh_triangle_area");
    
    geomlib_vec3 normal = geomlib_trimesh_triangle_normal(&mesh, 0);
    test(float_eq(normal.z, 1.0f), "trimesh_triangle_normal");
    
    geomlib_trimesh_free(&mesh);
}

static void test_trimesh_normals() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 64, 64);
    
    // Create a simple triangle
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 1.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    
    int result = geomlib_trimesh_compute_normals(&mesh);
    test(result == GEOMLIB_MESH_SUCCESS, "trimesh_compute_normals");
    test(mesh.normals != nullptr, "trimesh_normals_allocated");
    
    // All normals should point in Z direction for this flat triangle
    for (size_t i = 0; i < mesh.vertex_count; i++) {
        test(float_eq(mesh.normals[i].z, 1.0f), "trimesh_normal_direction");
    }
    
    geomlib_trimesh_free(&mesh);
}

static void test_trimesh_aabb() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 64, 64);
    
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(-1.0f, -2.0f, -3.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1.0f, 2.0f, 3.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 0.0f, 0.0f), nullptr, nullptr);
    
    geomlib_aabb box = geomlib_trimesh_compute_aabb(&mesh);
    
    test(float_eq(box.min.x, -1.0f) && float_eq(box.min.y, -2.0f) && float_eq(box.min.z, -3.0f), 
         "trimesh_aabb_min");
    test(float_eq(box.max.x, 1.0f) && float_eq(box.max.y, 2.0f) && float_eq(box.max.z, 3.0f), 
         "trimesh_aabb_max");
    
    geomlib_trimesh_free(&mesh);
}

static void test_aabb() {
    geomlib_aabb a = geomlib_aabb_create(
        geomlib_vec3_create(0.0f, 0.0f, 0.0f),
        geomlib_vec3_create(1.0f, 1.0f, 1.0f)
    );
    
    test(geomlib_aabb_contains_point(a, geomlib_vec3_create(0.5f, 0.5f, 0.5f)), 
         "aabb_contains_point_inside");
    test(!geomlib_aabb_contains_point(a, geomlib_vec3_create(2.0f, 0.5f, 0.5f)), 
         "aabb_contains_point_outside");
    
    geomlib_vec3 center = geomlib_aabb_center(a);
    test(float_eq(center.x, 0.5f) && float_eq(center.y, 0.5f) && float_eq(center.z, 0.5f), 
         "aabb_center");
    
    geomlib_vec3 size = geomlib_aabb_size(a);
    test(float_eq(size.x, 1.0f) && float_eq(size.y, 1.0f) && float_eq(size.z, 1.0f), 
         "aabb_size");
    
    float volume = geomlib_aabb_volume(a);
    test(float_eq(volume, 1.0f), "aabb_volume");
    
    float surface_area = geomlib_aabb_surface_area(a);
    test(float_eq(surface_area, 6.0f), "aabb_surface_area");
    
    // Test intersection
    geomlib_aabb b = geomlib_aabb_create(
        geomlib_vec3_create(0.5f, 0.5f, 0.5f),
        geomlib_vec3_create(1.5f, 1.5f, 1.5f)
    );
    test(geomlib_aabb_intersects(a, b), "aabb_intersects_true");
    
    geomlib_aabb c = geomlib_aabb_create(
        geomlib_vec3_create(2.0f, 2.0f, 2.0f),
        geomlib_vec3_create(3.0f, 3.0f, 3.0f)
    );
    test(!geomlib_aabb_intersects(a, c), "aabb_intersects_false");
}

static void test_create_box() {
    geomlib_trimesh mesh;
    int result = geomlib_create_box(&mesh, 2.0f, 2.0f, 2.0f);
    
    test(result == GEOMLIB_MESH_SUCCESS, "create_box_success");
    test(mesh.vertex_count == 24, "create_box_vertex_count");  // 4 per face * 6 faces
    test(mesh.triangle_count == 12, "create_box_triangle_count");  // 2 per face * 6 faces
    
    // Check surface area (6 faces * 4 = 24 for 2x2x2 box)
    float area = geomlib_trimesh_surface_area(&mesh);
    test(float_eq(area, 24.0f, 0.01f), "create_box_surface_area");
    
    geomlib_trimesh_free(&mesh);
}

static void test_create_sphere() {
    geomlib_trimesh mesh;
    int result = geomlib_create_sphere(&mesh, 1.0f, 16, 8);
    
    test(result == GEOMLIB_MESH_SUCCESS, "create_sphere_success");
    test(mesh.vertex_count > 0, "create_sphere_has_vertices");
    test(mesh.triangle_count > 0, "create_sphere_has_triangles");
    
    // Surface area of unit sphere should be approximately 4*pi
    float area = geomlib_trimesh_surface_area(&mesh);
    float expected_area = 4.0f * 3.14159265f;
    test(float_eq(area, expected_area, 0.5f), "create_sphere_surface_area");
    
    geomlib_trimesh_free(&mesh);
}

static void test_create_cylinder() {
    geomlib_trimesh mesh;
    int result = geomlib_create_cylinder(&mesh, 1.0f, 2.0f, 16);
    
    test(result == GEOMLIB_MESH_SUCCESS, "create_cylinder_success");
    test(mesh.vertex_count > 0, "create_cylinder_has_vertices");
    test(mesh.triangle_count > 0, "create_cylinder_has_triangles");
    
    geomlib_trimesh_free(&mesh);
}

static void test_create_plane() {
    geomlib_trimesh mesh;
    int result = geomlib_create_plane(&mesh, 2.0f, 2.0f, 2);
    
    test(result == GEOMLIB_MESH_SUCCESS, "create_plane_success");
    test(mesh.vertex_count == 9, "create_plane_vertex_count");  // 3x3 grid
    test(mesh.triangle_count == 8, "create_plane_triangle_count");  // 2x2 quads * 2
    
    float area = geomlib_trimesh_surface_area(&mesh);
    test(float_eq(area, 4.0f, 0.01f), "create_plane_surface_area");
    
    geomlib_trimesh_free(&mesh);
}

static void test_polymesh_triangulate() {
    geomlib_polymesh poly;
    geomlib_polymesh_init(&poly, 8, 4);
    
    // Create a simple quad
    geomlib_polymesh_add_vertex(&poly, geomlib_vec3_create(0.0f, 0.0f, 0.0f));
    geomlib_polymesh_add_vertex(&poly, geomlib_vec3_create(1.0f, 0.0f, 0.0f));
    geomlib_polymesh_add_vertex(&poly, geomlib_vec3_create(1.0f, 1.0f, 0.0f));
    geomlib_polymesh_add_vertex(&poly, geomlib_vec3_create(0.0f, 1.0f, 0.0f));
    
    uint32_t indices[] = {0, 1, 2, 3};
    geomlib_polymesh_add_face(&poly, indices, 4);
    
    geomlib_trimesh tri;
    int result = geomlib_polymesh_triangulate(&poly, &tri);
    
    test(result == GEOMLIB_MESH_SUCCESS, "polymesh_triangulate_success");
    test(tri.vertex_count == 4, "polymesh_triangulate_vertex_count");
    test(tri.triangle_count == 2, "polymesh_triangulate_triangle_count");
    
    geomlib_trimesh_free(&tri);
    geomlib_polymesh_free(&poly);
}

int main() {
    std::printf("=== Mesh Tests ===\n\n");
    
    test_polygon2d();
    test_polygon3d();
    test_trimesh_init();
    test_trimesh_triangle();
    test_trimesh_normals();
    test_trimesh_aabb();
    test_aabb();
    test_create_box();
    test_create_sphere();
    test_create_cylinder();
    test_create_plane();
    test_polymesh_triangulate();
    
    std::printf("\n=== Results: %d/%d tests passed ===\n", pass_count, test_count);
    
    return (pass_count == test_count) ? 0 : 1;
}
