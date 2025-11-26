/*
 * geomlib - Extended BVH Tests (Simplified)
 */

#include <cstdio>
#include <cmath>
#include "geomlib/geomlib_bvh.h"
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

static void test_single_triangle_bvh() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 3, 1);

    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1.0f, 0.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.0f, 1.0f, 0.0f), nullptr, nullptr);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);

    geomlib_bvh bvh;
    int result = geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);
    test(result == GEOMLIB_MESH_SUCCESS, "single_triangle_bvh_build");

    // Ray through center should hit
    geomlib_ray ray;
    ray.origin = geomlib_vec3_create(0.25f, 0.25f, 1.0f);
    ray.direction = geomlib_vec3_create(0.0f, 0.0f, -1.0f);
    ray.t_min = 0.0f;
    ray.t_max = 10.0f;

    geomlib_ray_hit hit;
    bool intersects = geomlib_bvh_trace_ray(&bvh, &ray, &hit);
    test(intersects && hit.hit, "single_triangle_bvh_ray_hit");

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

static void test_ray_miss() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 16, 8);

    geomlib_bvh bvh;
    geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);

    geomlib_ray ray;
    ray.origin = geomlib_vec3_create(0.0f, 0.0f, 5.0f);
    ray.direction = geomlib_vec3_create(0.0f, 1.0f, 0.0f);
    ray.t_min = 0.0f;
    ray.t_max = 100.0f;

    geomlib_ray_hit hit;
    bool intersects = geomlib_bvh_trace_ray(&bvh, &ray, &hit);
    test(!intersects || !hit.hit, "ray_miss_no_intersection");

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

static void test_bvh_quality_levels() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 16, 8);

    geomlib_bvh bvh1, bvh2, bvh3, bvh4;

    int r1 = geomlib_bvh_build(&bvh1, &mesh, GEOMLIB_BVH_QUALITY_FAST);
    int r2 = geomlib_bvh_build(&bvh2, &mesh, GEOMLIB_BVH_QUALITY_LOW);
    int r3 = geomlib_bvh_build(&bvh3, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);
    int r4 = geomlib_bvh_build(&bvh4, &mesh, GEOMLIB_BVH_QUALITY_HIGH);

    test(r1 == GEOMLIB_MESH_SUCCESS && r2 == GEOMLIB_MESH_SUCCESS &&
         r3 == GEOMLIB_MESH_SUCCESS && r4 == GEOMLIB_MESH_SUCCESS,
         "bvh_all_quality_levels_build");

    geomlib_bvh_free(&bvh1);
    geomlib_bvh_free(&bvh2);
    geomlib_bvh_free(&bvh3);
    geomlib_bvh_free(&bvh4);
    geomlib_trimesh_free(&mesh);
}

static void test_bvh_statistics() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 16, 8);

    geomlib_bvh bvh;
    geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);

    uint32_t max_depth, leaf_count;
    float avg_leaf_size;
    geomlib_bvh_get_stats(&bvh, &max_depth, &leaf_count, &avg_leaf_size);

    test(leaf_count > 0, "bvh_stats_has_leaves");
    test(max_depth >= 0, "bvh_stats_has_depth");
    test(avg_leaf_size > 0.0f, "bvh_stats_avg_leaf_size");

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

static void test_large_mesh_bvh() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 64, 32);

    test(mesh.triangle_count > 2000, "large_mesh_created");

    geomlib_bvh bvh;
    int result = geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);
    test(result == GEOMLIB_MESH_SUCCESS, "large_mesh_bvh_build");

    uint32_t max_depth, leaf_count;
    float avg_leaf_size;
    geomlib_bvh_get_stats(&bvh, &max_depth, &leaf_count, &avg_leaf_size);
    test(max_depth > 5, "large_mesh_bvh_has_depth");

    geomlib_ray ray;
    ray.origin = geomlib_vec3_create(0.0f, 0.0f, 3.0f);
    ray.direction = geomlib_vec3_create(0.0f, 0.0f, -1.0f);
    ray.t_min = 0.0f;
    ray.t_max = 10.0f;

    geomlib_ray_hit hit;
    bool intersects = geomlib_bvh_trace_ray(&bvh, &ray, &hit);
    test(intersects && hit.hit, "large_mesh_bvh_ray_intersection");

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

int main() {
    std::printf("=== Extended BVH Tests ===\n\n");

    test_single_triangle_bvh();
    test_ray_miss();
    test_bvh_quality_levels();
    test_bvh_statistics();
    test_large_mesh_bvh();

    std::printf("\n=== Results: %d/%d tests passed ===\n", pass_count, test_count);

    return (pass_count == test_count) ? 0 : 1;
}
