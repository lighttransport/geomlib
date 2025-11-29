/*
 * Comprehensive test for advanced geometry processing features
 */

#include <cstdio>
#include <cmath>
#include <cstring>
#include "geomlib/geomlib.h"

#define TEST_ASSERT(cond, msg) \
    if (!(cond)) { \
        printf("FAIL: %s (line %d)\n", msg, __LINE__); \
        return 1; \
    }

/* Test advancing front hole filling */
static int test_advancing_front() {
    printf("Testing advancing front triangulation...\n");
    
    /* Create a simple mesh with a hole */
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 10, 10);
    
    /* Create a square with a triangular hole */
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0, 0, 0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1, 0, 0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1, 1, 0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0, 1, 0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.5f, 0.5f, 0), NULL, NULL);
    
    /* Add triangles around the hole */
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 4);
    geomlib_trimesh_add_triangle(&mesh, 2, 3, 4);
    
    /* Detect hole */
    geomlib_holes holes = geomlib_detect_holes(&mesh);
    TEST_ASSERT(holes.hole_count > 0, "Should detect hole");
    
    /* Fill with advancing front */
    geomlib_fill_params params = geomlib_fill_params_default(GEOMLIB_FILL_ADVANCING_FRONT);
    int result = geomlib_fill_hole(&mesh, &holes.holes[0], &params);
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "Hole filling should succeed");
    
    geomlib_holes_free(&holes);
    geomlib_trimesh_free(&mesh);
    
    printf("  PASS\n");
    return 0;
}

/* Test Marching Cubes */
static int test_marching_cubes() {
    printf("Testing Marching Cubes...\n");

    /* Create a simple mesh and convert to SDF */
    geomlib_trimesh temp_mesh;
    geomlib_trimesh_init(&temp_mesh, 0, 0);
    geomlib_create_sphere(&temp_mesh, 0.5f, 8, 8);

    geomlib_sdf sdf = geomlib_mesh_to_sdf_resolution(&temp_mesh, 16);
    TEST_ASSERT(sdf.values != NULL, "SDF creation should succeed");

    /* Extract mesh with Marching Cubes */
    geomlib_trimesh mesh;
    int result = geomlib_sdf_to_mesh(&sdf, 0.0f, &mesh);
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "Marching Cubes should succeed");

    printf("  Generated %zu vertices, %zu triangles\n", mesh.vertex_count, mesh.triangle_count);

    geomlib_trimesh_free(&temp_mesh);
    geomlib_trimesh_free(&mesh);
    geomlib_sdf_free(&sdf);
    
    printf("  PASS\n");
    return 0;
}

/* Test BVH-accelerated self-intersection detection */
static int test_self_intersection() {
    printf("Testing BVH-accelerated self-intersection detection...\n");
    
    /* Create a mesh with self-intersection */
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 8, 4);
    
    /* First triangle */
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0, 0, 0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1, 0, 0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.5f, 1, 0), NULL, NULL);
    
    /* Second triangle (intersecting) */
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.5f, 0.3f, -0.1f), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.5f, 0.3f, 0.1f), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0.5f, 0.7f, 0), NULL, NULL);
    
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 3, 4, 5);
    
    /* Detect intersections */
    geomlib_self_intersection intersections = geomlib_detect_self_intersections(&mesh);
    
    printf("  Found %zu intersection pairs\n", intersections.pair_count);
    TEST_ASSERT(intersections.has_intersections, "Should detect intersections");
    
    geomlib_self_intersection_free(&intersections);
    geomlib_trimesh_free(&mesh);
    
    printf("  PASS\n");
    return 0;
}

/* Test nonmanifold vertex splitting */
static int test_nonmanifold_splitting() {
    printf("Testing nonmanifold vertex splitting...\n");
    
    /* Create a simple nonmanifold vertex (butterfly configuration) */
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 8, 4);
    
    /* Central vertex will be nonmanifold */
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0, 0, 0), NULL, NULL);  // 0: center
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1, 0, 0), NULL, NULL);  // 1
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0, 1, 0), NULL, NULL);  // 2
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(-1, 0, 0), NULL, NULL); // 3
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0, -1, 0), NULL, NULL); // 4
    
    /* Two separate triangle fans sharing only vertex 0 */
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 0, 3, 4);
    
    size_t orig_vertex_count = mesh.vertex_count;
    
    /* Split nonmanifold vertices */
    int splits = geomlib_split_nonmanifold_vertices(&mesh);
    
    printf("  Split %d nonmanifold vertices\n", splits);
    printf("  Vertices before: %zu, after: %zu\n", orig_vertex_count, mesh.vertex_count);
    
    geomlib_trimesh_free(&mesh);
    
    printf("  PASS\n");
    return 0;
}

/* Test surface reconstruction methods */
static int test_surface_reconstruction() {
    printf("Testing surface reconstruction (ball pivoting & alpha shape)...\n");
    
    /* Create a simple point cloud (cube corners) */
    geomlib_pointcloud pc;
    pc.point_count = 8;
    pc.points = (geomlib_vec3*)malloc(8 * sizeof(geomlib_vec3));
    pc.normals = NULL;
    
    /* Cube corners */
    pc.points[0] = geomlib_vec3_create(0, 0, 0);
    pc.points[1] = geomlib_vec3_create(1, 0, 0);
    pc.points[2] = geomlib_vec3_create(0, 1, 0);
    pc.points[3] = geomlib_vec3_create(1, 1, 0);
    pc.points[4] = geomlib_vec3_create(0, 0, 1);
    pc.points[5] = geomlib_vec3_create(1, 0, 1);
    pc.points[6] = geomlib_vec3_create(0, 1, 1);
    pc.points[7] = geomlib_vec3_create(1, 1, 1);
    
    /* Test ball pivoting */
    geomlib_surface_params params = geomlib_surface_params_default(GEOMLIB_SURFACE_BALL_PIVOTING);
    params.ball_radius = 1.5f;  /* Large enough to connect corners */
    
    geomlib_trimesh mesh;
    int result = geomlib_reconstruct_surface(&pc, &params, &mesh);
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "Ball pivoting should succeed");
    
    printf("  Ball pivoting: %zu triangles\n", mesh.triangle_count);
    geomlib_trimesh_free(&mesh);
    
    /* Test alpha shape */
    params = geomlib_surface_params_default(GEOMLIB_SURFACE_ALPHA_SHAPE);
    params.alpha = 1.5f;
    
    result = geomlib_reconstruct_surface(&pc, &params, &mesh);
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "Alpha shape should succeed");
    
    printf("  Alpha shape: %zu triangles\n", mesh.triangle_count);
    geomlib_trimesh_free(&mesh);
    
    free(pc.points);
    
    printf("  PASS\n");
    return 0;
}

/* Test SDF operations */
static int test_sdf_operations() {
    printf("Testing SDF operations (smoothing, blending, I/O)...\n");

    /* Create a simple mesh to convert to SDF */
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 0, 0);
    geomlib_create_box(&mesh, 1.0f, 1.0f, 1.0f);

    /* Create SDF from mesh */
    geomlib_sdf sdf1 = geomlib_mesh_to_sdf_resolution(&mesh, 16);
    TEST_ASSERT(sdf1.values != NULL, "SDF creation should succeed");

    geomlib_sdf sdf2 = geomlib_mesh_to_sdf_resolution(&mesh, 16);

    /* Test smoothing */
    int result = geomlib_sdf_smooth(&sdf1, 1.0f);
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "SDF smoothing should succeed");

    /* Test blending */
    geomlib_sdf blended = geomlib_sdf_blend(&sdf1, &sdf2, 0.5f);
    TEST_ASSERT(blended.values != NULL, "SDF blending should succeed");

    /* Test file I/O */
    result = geomlib_sdf_save(&sdf1, "/tmp/test_sdf.bin");
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "SDF save should succeed");

    geomlib_sdf loaded;
    result = geomlib_sdf_load(&loaded, "/tmp/test_sdf.bin");
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "SDF load should succeed");
    TEST_ASSERT(loaded.resolution[0] == sdf1.resolution[0], "Resolution should match");

    geomlib_trimesh_free(&mesh);
    geomlib_sdf_free(&sdf1);
    geomlib_sdf_free(&sdf2);
    geomlib_sdf_free(&blended);
    geomlib_sdf_free(&loaded);

    printf("  PASS\n");
    return 0;
}

int main() {
    printf("=== Advanced Features Test Suite ===\n\n");

    int failures = 0;

    failures += test_advancing_front();
    failures += test_marching_cubes();
    failures += test_self_intersection();
    failures += test_nonmanifold_splitting();
    failures += test_surface_reconstruction();
    failures += test_sdf_operations();
    
    printf("\n");
    if (failures == 0) {
        printf("=== ALL TESTS PASSED ===\n");
    } else {
        printf("=== %d TESTS FAILED ===\n", failures);
    }
    
    return failures;
}
