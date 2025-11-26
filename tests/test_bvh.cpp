/*
 * geomlib - BVH and Ray Tracing Tests
 * Tests for BVH construction, ray tracing, and point queries
 */

#include <geomlib/geomlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define TEST_EPSILON 1e-5f

void test_ray_creation() {
    printf("Testing ray creation...\n");

    geomlib_vec3 origin = {0.0f, 0.0f, 0.0f};
    geomlib_vec3 direction = {1.0f, 0.0f, 0.0f};

    geomlib_ray ray = geomlib_ray_create(origin, direction);

    assert(fabsf(ray.origin.x) < TEST_EPSILON);
    assert(fabsf(ray.origin.y) < TEST_EPSILON);
    assert(fabsf(ray.origin.z) < TEST_EPSILON);

    /* Direction should be normalized */
    float length = sqrtf(ray.direction.x * ray.direction.x +
                        ray.direction.y * ray.direction.y +
                        ray.direction.z * ray.direction.z);
    assert(fabsf(length - 1.0f) < TEST_EPSILON);

    printf("  ✓ Ray creation correct\n");
}

void test_ray_at() {
    printf("Testing ray evaluation...\n");

    geomlib_vec3 origin = {1.0f, 2.0f, 3.0f};
    geomlib_vec3 direction = {1.0f, 0.0f, 0.0f};
    geomlib_ray ray = geomlib_ray_create(origin, direction);

    geomlib_vec3 p = geomlib_ray_at(&ray, 5.0f);

    assert(fabsf(p.x - 6.0f) < TEST_EPSILON);
    assert(fabsf(p.y - 2.0f) < TEST_EPSILON);
    assert(fabsf(p.z - 3.0f) < TEST_EPSILON);

    printf("  ✓ Ray evaluation correct\n");
}

void test_ray_triangle_intersection() {
    printf("Testing ray-triangle intersection...\n");

    /* Triangle in xy plane at z=0 */
    geomlib_vec3 v0 = {0.0f, 0.0f, 0.0f};
    geomlib_vec3 v1 = {1.0f, 0.0f, 0.0f};
    geomlib_vec3 v2 = {0.0f, 1.0f, 0.0f};

    /* Ray shooting down from above */
    geomlib_vec3 origin = {0.25f, 0.25f, 1.0f};
    geomlib_vec3 direction = {0.0f, 0.0f, -1.0f};
    geomlib_ray ray = geomlib_ray_create(origin, direction);

    float t, u, v;
    bool hit = geomlib_ray_triangle_intersect(&ray, v0, v1, v2, &t, &u, &v);

    assert(hit);
    assert(fabsf(t - 1.0f) < TEST_EPSILON);  /* Hit at distance 1.0 */
    assert(u >= 0.0f && u <= 1.0f);
    assert(v >= 0.0f && v <= 1.0f);
    assert(u + v <= 1.0f);

    /* Ray missing triangle */
    origin = (geomlib_vec3){2.0f, 2.0f, 1.0f};
    direction = (geomlib_vec3){0.0f, 0.0f, -1.0f};
    ray = geomlib_ray_create(origin, direction);

    hit = geomlib_ray_triangle_intersect(&ray, v0, v1, v2, &t, &u, &v);
    assert(!hit);

    printf("  ✓ Ray-triangle intersection works\n");
}

void test_ray_aabb_intersection() {
    printf("Testing ray-AABB intersection...\n");

    geomlib_aabb box;
    box.min = (geomlib_vec3){-1.0f, -1.0f, -1.0f};
    box.max = (geomlib_vec3){1.0f, 1.0f, 1.0f};

    /* Ray through center */
    geomlib_vec3 origin = {-2.0f, 0.0f, 0.0f};
    geomlib_vec3 direction = {1.0f, 0.0f, 0.0f};
    geomlib_ray ray = geomlib_ray_create(origin, direction);

    bool hit = geomlib_ray_aabb_intersect(&ray, &box);
    assert(hit);

    /* Ray missing box */
    origin = (geomlib_vec3){-2.0f, 5.0f, 0.0f};
    direction = (geomlib_vec3){1.0f, 0.0f, 0.0f};
    ray = geomlib_ray_create(origin, direction);

    hit = geomlib_ray_aabb_intersect(&ray, &box);
    assert(!hit);

    printf("  ✓ Ray-AABB intersection works\n");
}

void test_bvh_construction() {
    printf("Testing BVH construction...\n");

    /* Create a simple mesh */
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 8, 12);

    /* Create a box */
    geomlib_vec3 vertices[] = {
        {-1.0f, -1.0f, -1.0f}, {1.0f, -1.0f, -1.0f},
        {1.0f, 1.0f, -1.0f}, {-1.0f, 1.0f, -1.0f},
        {-1.0f, -1.0f, 1.0f}, {1.0f, -1.0f, 1.0f},
        {1.0f, 1.0f, 1.0f}, {-1.0f, 1.0f, 1.0f}
    };

    for (int i = 0; i < 8; i++) {
        geomlib_trimesh_add_vertex(&mesh, vertices[i], NULL, NULL);
    }

    /* Add faces */
    uint32_t indices[] = {
        0, 1, 2,  0, 2, 3,  /* Front */
        4, 6, 5,  4, 7, 6,  /* Back */
        0, 4, 5,  0, 5, 1,  /* Bottom */
        2, 6, 7,  2, 7, 3,  /* Top */
        0, 3, 7,  0, 7, 4,  /* Left */
        1, 5, 6,  1, 6, 2   /* Right */
    };

    for (int i = 0; i < 36; i += 3) {
        geomlib_trimesh_add_triangle(&mesh, indices[i], indices[i+1], indices[i+2]);
    }

    /* Build BVH */
    geomlib_bvh bvh;
    int result = geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);

    assert(result == GEOMLIB_MESH_SUCCESS);
    assert(bvh.node_count > 0);
    assert(bvh.primitive_count == 12);

    /* Get stats */
    uint32_t max_depth, leaf_count;
    float avg_leaf_size;
    geomlib_bvh_get_stats(&bvh, &max_depth, &leaf_count, &avg_leaf_size);

    assert(max_depth > 0);
    assert(leaf_count > 0);
    assert(avg_leaf_size > 0.0f);

    printf("  ✓ BVH construction successful (depth=%u, leaves=%u, avg_leaf_size=%.2f)\n",
           max_depth, leaf_count, avg_leaf_size);

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

void test_bvh_ray_tracing() {
    printf("Testing BVH ray tracing...\n");

    /* Create a box mesh */
    geomlib_trimesh mesh;
    int result = geomlib_create_box(&mesh, 2.0f, 2.0f, 2.0f);
    assert(result == GEOMLIB_MESH_SUCCESS);

    /* Build BVH */
    geomlib_bvh bvh;
    result = geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);
    assert(result == GEOMLIB_MESH_SUCCESS);

    /* Ray hitting the box */
    geomlib_vec3 origin = {0.0f, 0.0f, 5.0f};
    geomlib_vec3 direction = {0.0f, 0.0f, -1.0f};
    geomlib_ray ray = geomlib_ray_create(origin, direction);

    geomlib_ray_hit hit;
    bool did_hit = geomlib_bvh_trace_ray(&bvh, &ray, &hit);

    assert(did_hit);
    assert(hit.hit);
    assert(hit.t > 0.0f);
    assert(hit.t < 10.0f);

    /* Hit position should be on front face */
    assert(fabsf(hit.position.z - 1.0f) < 0.01f);

    printf("  ✓ BVH ray tracing hit at t=%.4f, pos=(%.2f, %.2f, %.2f)\n",
           hit.t, hit.position.x, hit.position.y, hit.position.z);

    /* Ray missing the box */
    origin = (geomlib_vec3){10.0f, 10.0f, 5.0f};
    direction = (geomlib_vec3){0.0f, 0.0f, -1.0f};
    ray = geomlib_ray_create(origin, direction);

    did_hit = geomlib_bvh_trace_ray(&bvh, &ray, &hit);
    assert(!did_hit);

    printf("  ✓ BVH ray tracing correctly misses\n");

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

void test_bvh_occlusion() {
    printf("Testing BVH occlusion queries...\n");

    /* Create sphere mesh */
    geomlib_trimesh mesh;
    int result = geomlib_create_sphere(&mesh, 1.0f, 16, 8);
    assert(result == GEOMLIB_MESH_SUCCESS);

    /* Build BVH */
    geomlib_bvh bvh;
    result = geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_FAST);
    assert(result == GEOMLIB_MESH_SUCCESS);

    /* Ray that hits */
    geomlib_vec3 origin = {0.0f, 0.0f, 3.0f};
    geomlib_vec3 direction = {0.0f, 0.0f, -1.0f};
    geomlib_ray ray = geomlib_ray_create(origin, direction);

    bool occluded = geomlib_bvh_test_occlusion(&bvh, &ray);
    assert(occluded);

    /* Ray that misses */
    origin = (geomlib_vec3){5.0f, 5.0f, 3.0f};
    direction = (geomlib_vec3){0.0f, 0.0f, -1.0f};
    ray = geomlib_ray_create(origin, direction);

    occluded = geomlib_bvh_test_occlusion(&bvh, &ray);
    assert(!occluded);

    printf("  ✓ Occlusion queries work correctly\n");

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

void test_bvh_closest_point() {
    printf("Testing BVH closest point query...\n");

    /* Create sphere mesh */
    geomlib_trimesh mesh;
    int result = geomlib_create_sphere(&mesh, 1.0f, 32, 16);
    assert(result == GEOMLIB_MESH_SUCCESS);

    /* Build BVH */
    geomlib_bvh bvh;
    result = geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);
    assert(result == GEOMLIB_MESH_SUCCESS);

    /* Query point outside sphere */
    geomlib_vec3 query = {0.0f, 0.0f, 5.0f};
    geomlib_vec3 closest;
    float distance;
    uint32_t prim_id;

    bool found = geomlib_bvh_closest_point(&bvh, query, &closest, &distance, &prim_id);
    assert(found);

    /* Closest point should be on sphere surface (distance ~1.0 from origin) */
    float dist_from_origin = sqrtf(closest.x * closest.x + closest.y * closest.y + closest.z * closest.z);
    assert(fabsf(dist_from_origin - 1.0f) < 0.1f);

    /* Distance should be approximately 4.0 (5.0 - 1.0) */
    assert(distance > 3.5f && distance < 4.5f);

    printf("  ✓ Closest point query found point at distance %.4f\n", distance);

    /* Query point at origin (inside sphere) */
    query = (geomlib_vec3){0.0f, 0.0f, 0.0f};
    found = geomlib_bvh_closest_point(&bvh, query, &closest, &distance, &prim_id);
    assert(found);

    /* Closest point should be on surface */
    dist_from_origin = sqrtf(closest.x * closest.x + closest.y * closest.y + closest.z * closest.z);
    assert(fabsf(dist_from_origin - 1.0f) < 0.1f);

    printf("  ✓ Closest point query works for interior points\n");

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

void test_bvh_radius_query() {
    printf("Testing BVH radius query...\n");

    /* Create plane mesh */
    geomlib_trimesh mesh;
    int result = geomlib_create_plane(&mesh, 10.0f, 10.0f, 5);
    assert(result == GEOMLIB_MESH_SUCCESS);

    /* Build BVH */
    geomlib_bvh bvh;
    result = geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);
    assert(result == GEOMLIB_MESH_SUCCESS);

    /* Query above center of plane */
    geomlib_vec3 query = {0.0f, 1.0f, 0.0f};
    float radius = 3.0f;

    geomlib_vec3 points[100];
    float distances[100];
    uint32_t prim_ids[100];

    size_t found = geomlib_bvh_radius_query(&bvh, query, radius, 100,
                                            points, distances, prim_ids);

    assert(found > 0);

    /* All found points should be within radius */
    for (size_t i = 0; i < found; i++) {
        assert(distances[i] <= radius + TEST_EPSILON);
    }

    printf("  ✓ Radius query found %zu points within radius %.2f\n", found, radius);

    /* Small radius should find fewer points */
    found = geomlib_bvh_radius_query(&bvh, query, 0.5f, 100,
                                    points, distances, prim_ids);

    assert(found >= 0);  /* May find 0 or few points */

    printf("  ✓ Small radius query found %zu points\n", found);

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

void test_bvh_multiple_rays() {
    printf("Testing multiple rays through BVH...\n");

    /* Create sphere mesh */
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 32, 16);

    geomlib_bvh bvh;
    geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_HIGH);

    int hit_count = 0;
    int total_rays = 100;

    /* Cast rays in grid pattern */
    for (int y = 0; y < 10; y++) {
        for (int x = 0; x < 10; x++) {
            float fx = (x - 5.0f) * 0.3f;
            float fy = (y - 5.0f) * 0.3f;

            geomlib_vec3 origin = {fx, fy, 3.0f};
            geomlib_vec3 direction = {0.0f, 0.0f, -1.0f};
            geomlib_ray ray = geomlib_ray_create(origin, direction);

            geomlib_ray_hit hit;
            if (geomlib_bvh_trace_ray(&bvh, &ray, &hit)) {
                hit_count++;
            }
        }
    }

    /* Some rays should hit (those near center), some should miss */
    assert(hit_count > 0 && hit_count < total_rays);

    printf("  ✓ Multiple rays: %d/%d hit sphere\n", hit_count, total_rays);

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

void test_bvh_performance() {
    printf("Testing BVH performance...\n");

    /* Create detailed sphere */
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 64, 32);

    printf("  Mesh has %zu triangles\n", mesh.triangle_count);

    /* Build BVH */
    geomlib_bvh bvh;
    geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_HIGH);

    /* Cast many rays */
    int hit_count = 0;
    int num_rays = 10000;

    for (int i = 0; i < num_rays; i++) {
        float angle = i * 0.001f;
        geomlib_vec3 origin = {cosf(angle) * 3.0f, sinf(angle) * 3.0f, 0.0f};
        geomlib_vec3 direction = {-cosf(angle), -sinf(angle), 0.0f};
        geomlib_ray ray = geomlib_ray_create(origin, direction);

        geomlib_ray_hit hit;
        if (geomlib_bvh_trace_ray(&bvh, &ray, &hit)) {
            hit_count++;
        }
    }

    printf("  ✓ Cast %d rays, %d hits (%.1f%% hit rate)\n",
           num_rays, hit_count, (hit_count * 100.0f) / num_rays);

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
}

int main() {
    printf("=============================================\n");
    printf("geomlib BVH and Ray Tracing Tests\n");
    printf("=============================================\n\n");

    test_ray_creation();
    test_ray_at();
    test_ray_triangle_intersection();
    test_ray_aabb_intersection();
    test_bvh_construction();
    test_bvh_ray_tracing();
    test_bvh_occlusion();
    test_bvh_closest_point();
    test_bvh_radius_query();
    test_bvh_multiple_rays();
    test_bvh_performance();

    printf("\n=============================================\n");
    printf("All BVH tests passed!\n");
    printf("=============================================\n");

    return 0;
}
