/*
 * geomlib - Collision Detection and Mesh Reduction Tests
 */

#include <cstdio>
#include <cmath>
#include "geomlib/geomlib.h"

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

// Test OBB computation
static void test_obb_computation() {
    geomlib_trimesh mesh;
    geomlib_create_box(&mesh, 2.0f, 4.0f, 6.0f);

    geomlib_obb obb;
    int result = geomlib_compute_obb(&mesh, &obb);
    test(result == GEOMLIB_MESH_SUCCESS, "obb_computation_succeeds");

    // OBB volume should be positive and reasonable
    // (PCA-based OBB may not perfectly match original box volume)
    float obb_vol = geomlib_obb_volume(&obb);
    test(obb_vol > 0.0f && obb_vol < 1000.0f, "obb_volume_reasonable");

    // Center should be at origin
    test(float_eq(obb.center.x, 0.0f, 0.1f) &&
         float_eq(obb.center.y, 0.0f, 0.1f) &&
         float_eq(obb.center.z, 0.0f, 0.1f), "obb_center_at_origin");

    // Test point containment
    geomlib_vec3 inside = geomlib_vec3_create(0.5f, 0.5f, 0.5f);
    test(geomlib_obb_contains_point(&obb, inside), "obb_contains_inside_point");

    geomlib_vec3 outside = geomlib_vec3_create(10.0f, 10.0f, 10.0f);
    test(!geomlib_obb_contains_point(&obb, outside), "obb_excludes_outside_point");

    geomlib_trimesh_free(&mesh);
}

// Test OBB for sphere
static void test_obb_sphere() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 16, 8);

    geomlib_obb obb;
    int result = geomlib_compute_obb(&mesh, &obb);
    test(result == GEOMLIB_MESH_SUCCESS, "obb_sphere_succeeds");

    // For sphere, OBB should be roughly cubic
    test(float_eq(obb.half_extents.x, obb.half_extents.y, 0.2f) &&
         float_eq(obb.half_extents.y, obb.half_extents.z, 0.2f), "obb_sphere_roughly_cubic");

    geomlib_trimesh_free(&mesh);
}

// Test convex hull generation
static void test_convex_hull() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 16, 8);

    geomlib_trimesh hull;
    int result = geomlib_compute_convex_hull(&mesh, &hull);
    test(result == GEOMLIB_MESH_SUCCESS, "convex_hull_generation_succeeds");

    test(hull.vertex_count > 0 && hull.triangle_count > 0, "convex_hull_has_geometry");

    // Hull should enclose original mesh
    geomlib_aabb mesh_bounds = geomlib_trimesh_compute_aabb(&mesh);
    geomlib_aabb hull_bounds = geomlib_trimesh_compute_aabb(&hull);

    test(hull_bounds.min.x <= mesh_bounds.min.x + 0.1f &&
         hull_bounds.max.x >= mesh_bounds.max.x - 0.1f, "convex_hull_encloses_x");
    test(hull_bounds.min.y <= mesh_bounds.min.y + 0.1f &&
         hull_bounds.max.y >= mesh_bounds.max.y - 0.1f, "convex_hull_encloses_y");

    geomlib_trimesh_free(&hull);
    geomlib_trimesh_free(&mesh);
}

// Test convexity check
static void test_convexity_check() {
    // Box is convex
    geomlib_trimesh box;
    geomlib_create_box(&box, 1.0f, 1.0f, 1.0f);
    test(geomlib_is_convex(&box, 0.01f), "box_is_convex");
    geomlib_trimesh_free(&box);
}

// Test mesh simplification
static void test_mesh_simplification() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 32, 16);

    size_t original_count = mesh.triangle_count;
    test(original_count > 100, "mesh_has_triangles");

    // Simplify to 50% of original
    int result = geomlib_simplify_mesh_ratio(&mesh, 0.5f, 0);
    test(result == GEOMLIB_MESH_SUCCESS, "mesh_simplification_succeeds");

    test(mesh.triangle_count < original_count, "mesh_triangle_count_reduced");
    test(mesh.triangle_count >= 4, "mesh_has_minimum_triangles");

    geomlib_trimesh_free(&mesh);
}

// Test mesh simplification with target count
static void test_mesh_simplification_target() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 32, 16);

    size_t target = 50;
    int result = geomlib_simplify_mesh(&mesh, target, 0);
    test(result == GEOMLIB_MESH_SUCCESS, "mesh_simplification_target_succeeds");

    test(mesh.triangle_count <= target + 10, "mesh_close_to_target_count");

    geomlib_trimesh_free(&mesh);
}

// Test collision proxy generation - convex hull
static void test_proxy_convex_hull() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 32, 16);

    geomlib_trimesh proxy;
    int result = geomlib_generate_collision_proxy(&mesh, GEOMLIB_PROXY_CONVEX_HULL, 0, &proxy);
    test(result == GEOMLIB_MESH_SUCCESS, "proxy_convex_hull_succeeds");

    test(proxy.triangle_count > 0, "proxy_has_triangles");

    geomlib_trimesh_free(&proxy);
    geomlib_trimesh_free(&mesh);
}

// Test collision proxy generation - simplified
static void test_proxy_simplified() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 32, 16);

    size_t target = 50;
    geomlib_trimesh proxy;
    int result = geomlib_generate_collision_proxy(&mesh, GEOMLIB_PROXY_SIMPLIFIED, target, &proxy);
    test(result == GEOMLIB_MESH_SUCCESS, "proxy_simplified_succeeds");

    test(proxy.triangle_count > 0 && proxy.triangle_count <= target + 10,
         "proxy_simplified_count_reasonable");

    geomlib_trimesh_free(&proxy);
    geomlib_trimesh_free(&mesh);
}

// Test collision proxy generation - OBB
static void test_proxy_obb() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 16, 8);

    geomlib_trimesh proxy;
    int result = geomlib_generate_collision_proxy(&mesh, GEOMLIB_PROXY_OBB, 0, &proxy);
    test(result == GEOMLIB_MESH_SUCCESS, "proxy_obb_succeeds");

    // OBB proxy should be a box (12 triangles)
    test(proxy.triangle_count == 12, "proxy_obb_is_box");

    geomlib_trimesh_free(&proxy);
    geomlib_trimesh_free(&mesh);
}

// Test collision proxy generation - AABB
static void test_proxy_aabb() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 16, 8);

    geomlib_trimesh proxy;
    int result = geomlib_generate_collision_proxy(&mesh, GEOMLIB_PROXY_AABB, 0, &proxy);
    test(result == GEOMLIB_MESH_SUCCESS, "proxy_aabb_succeeds");

    // AABB proxy should be a box (12 triangles)
    test(proxy.triangle_count == 12, "proxy_aabb_is_box");

    geomlib_trimesh_free(&proxy);
    geomlib_trimesh_free(&mesh);
}

// Test box proxy generation
static void test_box_proxy() {
    geomlib_trimesh mesh;
    geomlib_create_cylinder(&mesh, 1.0f, 2.0f, 16);

    // AABB box
    geomlib_trimesh aabb_proxy;
    int result = geomlib_generate_box_proxy(&mesh, 0, &aabb_proxy);
    test(result == GEOMLIB_MESH_SUCCESS, "box_proxy_aabb_succeeds");
    test(aabb_proxy.triangle_count == 12, "box_proxy_aabb_count");

    // OBB box
    geomlib_trimesh obb_proxy;
    result = geomlib_generate_box_proxy(&mesh, 1, &obb_proxy);
    test(result == GEOMLIB_MESH_SUCCESS, "box_proxy_obb_succeeds");
    test(obb_proxy.triangle_count == 12, "box_proxy_obb_count");

    geomlib_trimesh_free(&obb_proxy);
    geomlib_trimesh_free(&aabb_proxy);
    geomlib_trimesh_free(&mesh);
}

// Test mesh quality analysis
static void test_mesh_quality() {
    geomlib_trimesh mesh;
    geomlib_create_box(&mesh, 1.0f, 1.0f, 1.0f);

    geomlib_mesh_quality quality;
    int result = geomlib_analyze_mesh_quality(&mesh, &quality);
    test(result == GEOMLIB_MESH_SUCCESS, "mesh_quality_analysis_succeeds");

    test(quality.min_edge_length > 0.0f, "quality_min_edge_positive");
    test(quality.max_edge_length >= quality.min_edge_length, "quality_max_ge_min");
    test(quality.avg_edge_length > 0.0f, "quality_avg_edge_positive");
    test(quality.degenerate_count == 0, "quality_no_degenerates");

    geomlib_trimesh_free(&mesh);
}

// Test collision cost estimation
static void test_collision_cost() {
    geomlib_trimesh simple, complex;
    geomlib_create_box(&simple, 1.0f, 1.0f, 1.0f);
    geomlib_create_sphere(&complex, 1.0f, 64, 32);

    float simple_cost = geomlib_estimate_collision_cost(&simple);
    float complex_cost = geomlib_estimate_collision_cost(&complex);

    test(simple_cost > 0.0f, "collision_cost_positive_simple");
    test(complex_cost > simple_cost, "collision_cost_complex_higher");

    geomlib_trimesh_free(&complex);
    geomlib_trimesh_free(&simple);
}

// Test proxy reduction comparison
static void test_proxy_comparison() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 64, 32);

    geomlib_trimesh aabb_proxy, simplified_proxy;

    geomlib_generate_collision_proxy(&mesh, GEOMLIB_PROXY_AABB, 0, &aabb_proxy);
    geomlib_generate_collision_proxy(&mesh, GEOMLIB_PROXY_SIMPLIFIED, 50, &simplified_proxy);

    float aabb_cost = geomlib_estimate_collision_cost(&aabb_proxy);
    float simplified_cost = geomlib_estimate_collision_cost(&simplified_proxy);
    float original_cost = geomlib_estimate_collision_cost(&mesh);

    test(aabb_cost < original_cost, "aabb_proxy_cheaper_than_original");
    test(simplified_cost < original_cost, "simplified_proxy_cheaper_than_original");

    geomlib_trimesh_free(&simplified_proxy);
    geomlib_trimesh_free(&aabb_proxy);
    geomlib_trimesh_free(&mesh);
}

// Test edge case: empty mesh
static void test_edge_cases() {
    geomlib_trimesh empty;
    geomlib_trimesh_init(&empty, 0, 0);

    geomlib_obb obb;
    test(geomlib_compute_obb(&empty, &obb) != GEOMLIB_MESH_SUCCESS,
         "obb_fails_on_empty_mesh");

    geomlib_trimesh hull;
    test(geomlib_compute_convex_hull(&empty, &hull) != GEOMLIB_MESH_SUCCESS,
         "convex_hull_fails_on_empty_mesh");

    geomlib_trimesh_free(&empty);
}

int main() {
    std::printf("=== Collision Detection and Mesh Reduction Tests ===\n\n");

    test_obb_computation();
    test_obb_sphere();
    test_convex_hull();
    test_convexity_check();
    test_mesh_simplification();
    test_mesh_simplification_target();
    test_proxy_convex_hull();
    test_proxy_simplified();
    test_proxy_obb();
    test_proxy_aabb();
    test_box_proxy();
    test_mesh_quality();
    test_collision_cost();
    test_proxy_comparison();
    test_edge_cases();

    std::printf("\n=== Results: %d/%d tests passed ===\n", pass_count, test_count);

    return (pass_count == test_count) ? 0 : 1;
}
