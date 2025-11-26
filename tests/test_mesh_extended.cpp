/*
 * geomlib - Extended Mesh Tests
 * Tests edge cases, volume/centroid, and stress scenarios
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

// Test degenerate triangles
static void test_degenerate_triangles() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 64, 64);

    // Degenerate triangle (all points collinear)
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(2.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);

    float area = geomlib_trimesh_triangle_area(&mesh, 0);
    test(float_eq(area, 0.0f, 1e-6f), "degenerate_triangle_zero_area");

    // Zero-length edge triangle
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_triangle(&mesh, 3, 4, 5);

    area = geomlib_trimesh_triangle_area(&mesh, 1);
    test(float_eq(area, 0.0f, 1e-6f), "zero_length_edge_triangle");

    geomlib_trimesh_free(&mesh);
}

// Test mesh volume calculation
static void test_mesh_volume() {
    geomlib_trimesh mesh;

    // Unit cube should have volume 1
    geomlib_create_box(&mesh, 1.0f, 1.0f, 1.0f);
    float volume = geomlib_trimesh_volume(&mesh);
    test(float_eq(volume, 1.0f, 0.01f), "unit_cube_volume");
    geomlib_trimesh_free(&mesh);

    // 2x2x2 cube should have volume 8
    geomlib_create_box(&mesh, 2.0f, 2.0f, 2.0f);
    volume = geomlib_trimesh_volume(&mesh);
    test(float_eq(volume, 8.0f, 0.01f), "2x2x2_cube_volume");
    geomlib_trimesh_free(&mesh);

    // Sphere volume: 4/3 * pi * r^3
    geomlib_create_sphere(&mesh, 1.0f, 32, 16);
    volume = geomlib_trimesh_volume(&mesh);
    float expected = (4.0f / 3.0f) * 3.14159265f;
    test(float_eq(volume, expected, 0.1f), "unit_sphere_volume");
    geomlib_trimesh_free(&mesh);
}

// Test mesh centroid calculation
static void test_mesh_centroid() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 64, 64);

    // Symmetric triangle centered at origin
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(-1.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 1.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);

    geomlib_vec3 centroid = geomlib_trimesh_centroid(&mesh);
    test(float_eq(centroid.x, 0.0f, 0.01f) && float_eq(centroid.y, 0.333f, 0.01f),
         "triangle_centroid");

    geomlib_trimesh_free(&mesh);

    // Box centroid should be at center
    geomlib_create_box(&mesh, 2.0f, 4.0f, 6.0f);
    centroid = geomlib_trimesh_centroid(&mesh);
    test(float_eq(centroid.x, 0.0f, 0.01f) &&
         float_eq(centroid.y, 0.0f, 0.01f) &&
         float_eq(centroid.z, 0.0f, 0.01f), "box_centroid");
    geomlib_trimesh_free(&mesh);
}

// Test reserve and capacity
static void test_reserve_capacity() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 8, 8);

    test(mesh.vertex_capacity >= 8, "initial_vertex_capacity");
    test(mesh.triangle_capacity >= 8, "initial_triangle_capacity");

    int result = geomlib_trimesh_reserve_vertices(&mesh, 100);
    test(result == GEOMLIB_MESH_SUCCESS && mesh.vertex_capacity >= 100,
         "reserve_vertices_expands");

    result = geomlib_trimesh_reserve_triangles(&mesh, 200);
    test(result == GEOMLIB_MESH_SUCCESS && mesh.triangle_capacity >= 200,
         "reserve_triangles_expands");

    geomlib_trimesh_free(&mesh);
}

// Test clear operations
static void test_clear_operations() {
    geomlib_trimesh mesh;
    geomlib_create_box(&mesh, 1.0f, 1.0f, 1.0f);

    test(mesh.vertex_count > 0 && mesh.triangle_count > 0, "mesh_has_data");

    geomlib_trimesh_clear(&mesh);
    test(mesh.vertex_count == 0 && mesh.triangle_count == 0, "mesh_cleared");

    // Should be able to reuse after clear
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 0.0f, 0.0f), nullptr, nullptr);
    test(mesh.vertex_count == 1, "mesh_reusable_after_clear");

    geomlib_trimesh_free(&mesh);
}

// Test triangle vertex retrieval
static void test_triangle_vertices() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 64, 64);

    geomlib_vec3 v0 = geomlib_vec3_create(1.0f, 2.0f, 3.0f);
    geomlib_vec3 v1 = geomlib_vec3_create(4.0f, 5.0f, 6.0f);
    geomlib_vec3 v2 = geomlib_vec3_create(7.0f, 8.0f, 9.0f);

    geomlib_trimesh_add_vertex(&mesh, v0, nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, v1, nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, v2, nullptr, nullptr);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);

    geomlib_vec3 out0, out1, out2;
    geomlib_trimesh_get_triangle_vertices(&mesh, 0, &out0, &out1, &out2);

    test(float_eq(out0.x, v0.x) && float_eq(out0.y, v0.y) && float_eq(out0.z, v0.z),
         "triangle_vertex_0");
    test(float_eq(out1.x, v1.x) && float_eq(out1.y, v1.y) && float_eq(out1.z, v1.z),
         "triangle_vertex_1");
    test(float_eq(out2.x, v2.x) && float_eq(out2.y, v2.y) && float_eq(out2.z, v2.z),
         "triangle_vertex_2");

    geomlib_trimesh_free(&mesh);
}

// Test extreme coordinates
static void test_extreme_coordinates() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 64, 64);

    // Very large coordinates
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1e6f, 1e6f, 1e6f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1e6f + 1.0f, 1e6f, 1e6f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1e6f, 1e6f + 1.0f, 1e6f), nullptr, nullptr);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);

    float area = geomlib_trimesh_triangle_area(&mesh, 0);
    test(float_eq(area, 0.5f, 0.1f), "large_coordinate_triangle_area");

    geomlib_trimesh_free(&mesh);

    // Very small coordinates
    geomlib_trimesh_init(&mesh, 64, 64);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1e-3f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 1e-3f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);

    area = geomlib_trimesh_triangle_area(&mesh, 0);
    test(float_eq(area, 0.5e-6f, 1e-7f), "small_coordinate_triangle_area");

    geomlib_trimesh_free(&mesh);
}

// Test AABB merge and expand
static void test_aabb_merge_expand() {
    geomlib_aabb a = geomlib_aabb_create(
        geomlib_vec3_create(0.0f, 0.0f, 0.0f),
        geomlib_vec3_create(1.0f, 1.0f, 1.0f)
    );

    geomlib_aabb b = geomlib_aabb_create(
        geomlib_vec3_create(0.5f, 0.5f, 0.5f),
        geomlib_vec3_create(2.0f, 2.0f, 2.0f)
    );

    geomlib_aabb merged = geomlib_aabb_merge(a, b);
    test(float_eq(merged.min.x, 0.0f) && float_eq(merged.max.x, 2.0f), "aabb_merge_x");
    test(float_eq(merged.min.y, 0.0f) && float_eq(merged.max.y, 2.0f), "aabb_merge_y");
    test(float_eq(merged.min.z, 0.0f) && float_eq(merged.max.z, 2.0f), "aabb_merge_z");

    // Test expand
    geomlib_aabb expanded = geomlib_aabb_expand(a, geomlib_vec3_create(5.0f, 5.0f, 5.0f));
    test(float_eq(expanded.max.x, 5.0f) && float_eq(expanded.max.y, 5.0f) &&
         float_eq(expanded.max.z, 5.0f), "aabb_expand_point");
}

// Test AABB closest point
static void test_aabb_closest_point() {
    geomlib_aabb box = geomlib_aabb_create(
        geomlib_vec3_create(0.0f, 0.0f, 0.0f),
        geomlib_vec3_create(1.0f, 1.0f, 1.0f)
    );

    // Point inside - should return the point itself
    geomlib_vec3 p1 = geomlib_vec3_create(0.5f, 0.5f, 0.5f);
    geomlib_vec3 closest = geomlib_aabb_closest_point(&box, p1);
    test(float_eq(closest.x, 0.5f) && float_eq(closest.y, 0.5f) && float_eq(closest.z, 0.5f),
         "aabb_closest_point_inside");

    // Point outside on X axis
    geomlib_vec3 p2 = geomlib_vec3_create(2.0f, 0.5f, 0.5f);
    closest = geomlib_aabb_closest_point(&box, p2);
    test(float_eq(closest.x, 1.0f) && float_eq(closest.y, 0.5f) && float_eq(closest.z, 0.5f),
         "aabb_closest_point_outside_x");

    // Point outside at corner
    geomlib_vec3 p3 = geomlib_vec3_create(2.0f, 2.0f, 2.0f);
    closest = geomlib_aabb_closest_point(&box, p3);
    test(float_eq(closest.x, 1.0f) && float_eq(closest.y, 1.0f) && float_eq(closest.z, 1.0f),
         "aabb_closest_point_corner");

    // Point below min
    geomlib_vec3 p4 = geomlib_vec3_create(-1.0f, -1.0f, -1.0f);
    closest = geomlib_aabb_closest_point(&box, p4);
    test(float_eq(closest.x, 0.0f) && float_eq(closest.y, 0.0f) && float_eq(closest.z, 0.0f),
         "aabb_closest_point_below_min");
}

// Test empty AABB
static void test_aabb_empty() {
    geomlib_aabb empty = geomlib_aabb_empty();
    test(empty.min.x > empty.max.x, "aabb_empty_inverted");

    // Expanding empty should give point
    geomlib_vec3 p = geomlib_vec3_create(1.0f, 2.0f, 3.0f);
    geomlib_aabb expanded = geomlib_aabb_expand(empty, p);
    test(float_eq(expanded.min.x, 1.0f) && float_eq(expanded.max.x, 1.0f) &&
         float_eq(expanded.min.y, 2.0f) && float_eq(expanded.max.y, 2.0f) &&
         float_eq(expanded.min.z, 3.0f) && float_eq(expanded.max.z, 3.0f),
         "aabb_empty_expand_to_point");
}

// Test polygon2d edge cases
static void test_polygon2d_edge_cases() {
    geomlib_polygon2d poly;
    geomlib_polygon2d_init(&poly, 4);

    // Triangle (3 vertices)
    geomlib_polygon2d_add_vertex(&poly, geomlib_vec2_create(0.0f, 0.0f));
    geomlib_polygon2d_add_vertex(&poly, geomlib_vec2_create(1.0f, 0.0f));
    geomlib_polygon2d_add_vertex(&poly, geomlib_vec2_create(0.0f, 1.0f));

    float area = geomlib_polygon2d_area(&poly);
    test(float_eq(area, 0.5f), "polygon2d_triangle_area");

    test(geomlib_polygon2d_is_convex(&poly), "polygon2d_triangle_is_convex");

    // Test boundary point
    test(geomlib_polygon2d_contains_point(&poly, geomlib_vec2_create(0.0f, 0.0f)),
         "polygon2d_contains_vertex");

    geomlib_polygon2d_free(&poly);
}

// Test large mesh (stress test)
static void test_large_mesh() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 64, 32);  // High resolution sphere

    test(mesh.vertex_count > 1000, "large_mesh_vertex_count");
    test(mesh.triangle_count > 1000, "large_mesh_triangle_count");

    // Should still compute properties correctly
    geomlib_aabb box = geomlib_trimesh_compute_aabb(&mesh);
    geomlib_vec3 center = geomlib_aabb_center(box);
    test(float_eq(center.x, 0.0f, 0.1f) && float_eq(center.y, 0.0f, 0.1f) &&
         float_eq(center.z, 0.0f, 0.1f), "large_mesh_aabb_centered");

    float surface_area = geomlib_trimesh_surface_area(&mesh);
    float expected = 4.0f * 3.14159265f;
    test(float_eq(surface_area, expected, 0.2f), "large_mesh_surface_area");

    geomlib_trimesh_free(&mesh);
}

// Test primitive generation with different parameters
static void test_primitive_parameters() {
    geomlib_trimesh mesh;

    // Low resolution sphere
    geomlib_create_sphere(&mesh, 2.0f, 8, 4);
    test(mesh.vertex_count > 0 && mesh.triangle_count > 0, "low_res_sphere");
    geomlib_trimesh_free(&mesh);

    // Flat cylinder
    geomlib_create_cylinder(&mesh, 5.0f, 0.1f, 32);
    test(mesh.vertex_count > 0 && mesh.triangle_count > 0, "flat_cylinder");
    geomlib_trimesh_free(&mesh);

    // Large plane with subdivisions
    geomlib_create_plane(&mesh, 100.0f, 100.0f, 10);
    test(mesh.vertex_count == 11 * 11, "plane_subdivision_vertices");
    test(mesh.triangle_count == 10 * 10 * 2, "plane_subdivision_triangles");
    geomlib_trimesh_free(&mesh);
}

int main() {
    std::printf("=== Extended Mesh Tests ===\n\n");

    test_degenerate_triangles();
    test_mesh_volume();
    test_mesh_centroid();
    test_reserve_capacity();
    test_clear_operations();
    test_triangle_vertices();
    test_extreme_coordinates();
    test_aabb_merge_expand();
    test_aabb_closest_point();
    test_aabb_empty();
    test_polygon2d_edge_cases();
    test_large_mesh();
    test_primitive_parameters();

    std::printf("\n=== Results: %d/%d tests passed ===\n", pass_count, test_count);

    return (pass_count == test_count) ? 0 : 1;
}
