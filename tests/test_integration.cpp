/*
 * geomlib - Integration Tests (Simplified)
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

static void test_mesh_montecarlo_integration() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 16, 8);

    geomlib_mesh_sample_distribution dist;
    int result = geomlib_mesh_build_sample_distribution(&mesh, &dist);
    test(result == GEOMLIB_MESH_SUCCESS, "mesh_montecarlo_build_distribution");

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 12345, 0);

    geomlib_vec3 point = geomlib_sample_mesh_surface(&mesh, &dist, &rng, nullptr);
    float dist_from_center = geomlib_vec3_length(point);
    test(float_eq(dist_from_center, 1.0f, 0.15f), "sampled_point_on_sphere_surface");

    geomlib_mesh_free_sample_distribution(&dist);
    geomlib_trimesh_free(&mesh);
}

static void test_raytrace_primitives_integration() {
    geomlib_trimesh box, sphere;
    geomlib_create_box(&box, 1.0f, 1.0f, 1.0f);
    geomlib_create_sphere(&sphere, 0.5f, 16, 8);

    geomlib_bvh bvh_box, bvh_sphere;
    geomlib_bvh_build(&bvh_box, &box, GEOMLIB_BVH_QUALITY_MEDIUM);
    geomlib_bvh_build(&bvh_sphere, &sphere, GEOMLIB_BVH_QUALITY_MEDIUM);

    geomlib_ray ray;
    ray.origin = geomlib_vec3_create(0.0f, 0.0f, 3.0f);
    ray.direction = geomlib_vec3_create(0.0f, 0.0f, -1.0f);
    ray.t_min = 0.0f;
    ray.t_max = 10.0f;

    geomlib_ray_hit hit;
    bool hit_box = geomlib_bvh_trace_ray(&bvh_box, &ray, &hit);
    test(hit_box && hit.hit, "raytrace_box_hit");

    bool hit_sphere = geomlib_bvh_trace_ray(&bvh_sphere, &ray, &hit);
    test(hit_sphere && hit.hit, "raytrace_sphere_hit");

    geomlib_bvh_free(&bvh_box);
    geomlib_bvh_free(&bvh_sphere);
    geomlib_trimesh_free(&box);
    geomlib_trimesh_free(&sphere);
}

static void test_random_montecarlo_consistency() {
    geomlib_pcg32 rng1, rng2;
    geomlib_pcg32_init(&rng1, 999, 0);
    geomlib_pcg32_init(&rng2, 999, 0);

    geomlib_vec3 normal = geomlib_vec3_create(0.0f, 0.0f, 1.0f);
    geomlib_vec3 s1 = geomlib_sample_hemisphere_cosine(normal, &rng1);
    geomlib_vec3 s2 = geomlib_sample_hemisphere_cosine(normal, &rng2);

    test(float_eq(s1.x, s2.x) && float_eq(s1.y, s2.y) && float_eq(s1.z, s2.z),
         "random_montecarlo_consistent_samples");
}

static void test_mesh_area_montecarlo_integration() {
    geomlib_trimesh mesh;
    geomlib_create_box(&mesh, 2.0f, 2.0f, 2.0f);

    float exact_area = geomlib_trimesh_surface_area(&mesh);
    test(float_eq(exact_area, 24.0f, 0.01f), "mesh_exact_surface_area");

    geomlib_mesh_sample_distribution dist;
    geomlib_mesh_build_sample_distribution(&mesh, &dist);

    test(float_eq(dist.total_area, exact_area, 0.01f), "montecarlo_area_matches_exact");

    geomlib_mesh_free_sample_distribution(&dist);
    geomlib_trimesh_free(&mesh);
}

static void test_linalg_mesh_transform() {
    geomlib_trimesh mesh;
    geomlib_create_box(&mesh, 1.0f, 1.0f, 1.0f);

    geomlib_mat4 transform = geomlib_mat4_scale_xyz(2.0f, 2.0f, 2.0f);

    for (size_t i = 0; i < mesh.vertex_count; i++) {
        geomlib_vec3 v = mesh.vertices[i];
        geomlib_vec4 v4 = geomlib_vec4_create(v.x, v.y, v.z, 1.0f);
        geomlib_vec4 transformed = geomlib_mat4_mul_vec4(transform, v4);
        mesh.vertices[i] = geomlib_vec3_create(transformed.x, transformed.y, transformed.z);
    }

    float new_area = geomlib_trimesh_surface_area(&mesh);
    test(float_eq(new_area, 24.0f, 0.1f), "linalg_mesh_scaled_area");

    float new_volume = geomlib_trimesh_volume(&mesh);
    test(float_eq(new_volume, 8.0f, 0.1f), "linalg_mesh_scaled_volume");

    geomlib_trimesh_free(&mesh);
}

static void test_aabb_mesh_bounds() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 16, 8);

    geomlib_aabb bounds = geomlib_trimesh_compute_aabb(&mesh);

    test(float_eq(bounds.min.x, -1.0f, 0.1f) && float_eq(bounds.max.x, 1.0f, 0.1f),
         "aabb_sphere_bounds_x");

    geomlib_vec3 centroid = geomlib_trimesh_centroid(&mesh);
    test(geomlib_aabb_contains_point(bounds, centroid), "aabb_contains_centroid");

    geomlib_trimesh_free(&mesh);
}

static void test_full_pipeline() {
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 32, 16);
    test(mesh.triangle_count > 0, "pipeline_mesh_created");

    geomlib_bvh bvh;
    int result = geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);
    test(result == GEOMLIB_MESH_SUCCESS, "pipeline_bvh_built");

    geomlib_mesh_sample_distribution dist;
    result = geomlib_mesh_build_sample_distribution(&mesh, &dist);
    test(result == GEOMLIB_MESH_SUCCESS, "pipeline_sampling_setup");

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    geomlib_vec3 point = geomlib_sample_mesh_surface(&mesh, &dist, &rng, nullptr);
    
    geomlib_ray ray;
    ray.origin = geomlib_vec3_create(0.0f, 0.0f, 0.0f);
    ray.direction = geomlib_vec3_normalize(point);
    ray.t_min = 0.0f;
    ray.t_max = 2.0f;

    geomlib_ray_hit hit;
    bool intersects = geomlib_bvh_trace_ray(&bvh, &ray, &hit);
    test(intersects, "pipeline_ray_hits");

    geomlib_mesh_free_sample_distribution(&dist);
    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

static void test_polygon_mesh_bvh() {
    geomlib_polymesh poly;
    geomlib_polymesh_init(&poly, 4, 1);

    geomlib_polymesh_add_vertex(&poly, geomlib_vec3_create(0.0f, 0.0f, 0.0f));
    geomlib_polymesh_add_vertex(&poly, geomlib_vec3_create(1.0f, 0.0f, 0.0f));
    geomlib_polymesh_add_vertex(&poly, geomlib_vec3_create(1.0f, 1.0f, 0.0f));
    geomlib_polymesh_add_vertex(&poly, geomlib_vec3_create(0.0f, 1.0f, 0.0f));

    uint32_t face[] = {0, 1, 2, 3};
    geomlib_polymesh_add_face(&poly, face, 4);

    geomlib_trimesh tri;
    int result = geomlib_polymesh_triangulate(&poly, &tri);
    test(result == GEOMLIB_MESH_SUCCESS, "polygon_triangulation");

    geomlib_bvh bvh;
    result = geomlib_bvh_build(&bvh, &tri, GEOMLIB_BVH_QUALITY_MEDIUM);
    test(result == GEOMLIB_MESH_SUCCESS, "polygon_mesh_bvh_build");

    geomlib_ray ray;
    ray.origin = geomlib_vec3_create(0.5f, 0.5f, 1.0f);
    ray.direction = geomlib_vec3_create(0.0f, 0.0f, -1.0f);
    ray.t_min = 0.0f;
    ray.t_max = 10.0f;

    geomlib_ray_hit hit;
    bool intersects = geomlib_bvh_trace_ray(&bvh, &ray, &hit);
    test(intersects && hit.hit, "polygon_mesh_bvh_ray_hit");

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&tri);
    geomlib_polymesh_free(&poly);
}

int main() {
    std::printf("=== Integration Tests ===\n\n");

    test_mesh_montecarlo_integration();
    test_raytrace_primitives_integration();
    test_random_montecarlo_consistency();
    test_mesh_area_montecarlo_integration();
    test_linalg_mesh_transform();
    test_aabb_mesh_bounds();
    test_full_pipeline();
    test_polygon_mesh_bvh();

    std::printf("\n=== Results: %d/%d tests passed ===\n", pass_count, test_count);

    return (pass_count == test_count) ? 0 : 1;
}
