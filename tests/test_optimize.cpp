/*
 * Tests for mesh optimization features
 */

#include <cstdio>
#include <cmath>
#include "geomlib/geomlib.h"

#define TEST_ASSERT(cond, msg) \
    if (!(cond)) { \
        printf("FAIL: %s (line %d)\n", msg, __LINE__); \
        return 1; \
    }

/* Test triangle strip generation */
static int test_triangle_strips() {
    printf("Testing triangle strip generation...\n");

    /* Create a simple mesh (two connected triangles forming a quad) */
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 4, 2);

    /* Add vertices for a quad */
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0, 0, 0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1, 0, 0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1, 1, 0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0, 1, 0), NULL, NULL);

    /* Add two triangles */
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 0, 2, 3);

    /* Generate strips */
    geomlib_tristrip_group strips = geomlib_generate_tristrips(&mesh, false);

    printf("  Generated %zu strip(s)\n", strips.strip_count);
    TEST_ASSERT(strips.strip_count > 0, "Should generate at least one strip");

    /* Verify strips cover all triangles */
    size_t total_triangles = 0;
    for (size_t i = 0; i < strips.strip_count; i++) {
        TEST_ASSERT(strips.strips[i].index_count >= 3, "Strip should have at least 3 indices");
        /* Each strip with N indices represents N-2 triangles */
        total_triangles += strips.strips[i].index_count - 2;
        printf("  Strip %zu: %zu indices (%zu triangles)\n",
               i, strips.strips[i].index_count, strips.strips[i].index_count - 2);
    }

    TEST_ASSERT(total_triangles == mesh.triangle_count, "Strips should cover all triangles");

    geomlib_tristrip_group_free(&strips);
    geomlib_trimesh_free(&mesh);

    printf("  PASS\n");
    return 0;
}

/* Test index compression */
static int test_index_compression() {
    printf("Testing index compression...\n");

    /* Create a simple index sequence */
    uint32_t indices[] = {0, 1, 2, 2, 1, 3, 3, 1, 4, 4, 1, 5};
    size_t index_count = sizeof(indices) / sizeof(indices[0]);

    /* Compress */
    geomlib_compressed_indices compressed = geomlib_compress_indices(
        indices, index_count, 6
    );

    TEST_ASSERT(compressed.data != NULL, "Compression should succeed");
    TEST_ASSERT(compressed.original_count == index_count, "Original count should match");

    /* Check compression ratio */
    size_t original_size = index_count * sizeof(uint32_t);
    printf("  Original size: %zu bytes\n", original_size);
    printf("  Compressed size: %zu bytes\n", compressed.compressed_size);
    printf("  Compression ratio: %.2f%%\n",
           100.0f * compressed.compressed_size / original_size);

    TEST_ASSERT(compressed.compressed_size <= original_size,
                "Compressed size should not exceed original");

    /* Decompress */
    uint32_t decompressed[12];
    int result = geomlib_decompress_indices(&compressed, decompressed);
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "Decompression should succeed");

    /* Verify data matches */
    for (size_t i = 0; i < index_count; i++) {
        TEST_ASSERT(decompressed[i] == indices[i], "Decompressed data should match");
    }

    geomlib_compressed_indices_free(&compressed);

    const uint32_t wide[] = {0, UINT32_MAX, 0x80000000u, 0, 400};
    compressed = geomlib_compress_indices(wide, 5, UINT32_MAX);
    TEST_ASSERT(geomlib_decompress_indices(&compressed, decompressed) == GEOMLIB_MESH_SUCCESS,
                "Full-width deltas should round-trip");
    for (size_t i = 0; i < 5; i++)
        TEST_ASSERT(decompressed[i] == wide[i], "Full-width delta must retain its bits");
    for (size_t size = 0; size < compressed.compressed_size; size++) {
        auto truncated = compressed;
        truncated.compressed_size = size;
        TEST_ASSERT(geomlib_decompress_indices(&truncated, decompressed) == GEOMLIB_MESH_ERROR_INVALID,
                    "Truncated streams must fail");
    }
    geomlib_compressed_indices_free(&compressed);

    printf("  PASS\n");
    return 0;
}

/* Test compression with larger mesh */
static int test_compression_sphere() {
    printf("Testing compression on sphere mesh...\n");

    /* Create a sphere */
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 16, 16);

    printf("  Sphere: %zu vertices, %zu triangles\n",
           mesh.vertex_count, mesh.triangle_count);

    /* Compress indices */
    geomlib_compressed_indices compressed = geomlib_compress_indices(
        mesh.indices, mesh.triangle_count * 3, (uint32_t)mesh.vertex_count
    );

    TEST_ASSERT(compressed.data != NULL, "Compression should succeed");

    size_t original_size = mesh.triangle_count * 3 * sizeof(uint32_t);
    printf("  Original size: %zu bytes\n", original_size);
    printf("  Compressed size: %zu bytes\n", compressed.compressed_size);
    printf("  Compression ratio: %.2f%%\n",
           100.0f * compressed.compressed_size / original_size);

    /* Decompress and verify */
    uint32_t *decompressed = (uint32_t*)malloc(mesh.triangle_count * 3 * sizeof(uint32_t));
    int result = geomlib_decompress_indices(&compressed, decompressed);
    TEST_ASSERT(result == GEOMLIB_MESH_SUCCESS, "Decompression should succeed");

    /* Verify first few indices */
    for (size_t i = 0; i < 10 && i < mesh.triangle_count * 3; i++) {
        TEST_ASSERT(decompressed[i] == mesh.indices[i],
                    "Decompressed indices should match");
    }

    free(decompressed);
    geomlib_compressed_indices_free(&compressed);
    geomlib_trimesh_free(&mesh);

    printf("  PASS\n");
    return 0;
}

/* Test strips on larger mesh */
static int test_strips_sphere() {
    printf("Testing triangle strips on sphere...\n");

    /* Create a sphere */
    geomlib_trimesh mesh;
    geomlib_create_sphere(&mesh, 1.0f, 8, 8);

    printf("  Sphere: %zu triangles\n", mesh.triangle_count);

    /* Generate strips */
    geomlib_tristrip_group strips = geomlib_generate_tristrips(&mesh, false);

    printf("  Generated %zu strip(s)\n", strips.strip_count);
    TEST_ASSERT(strips.strip_count > 0, "Should generate strips");

    /* Count total triangles in strips */
    size_t total_triangles = 0;
    size_t total_indices = 0;
    for (size_t i = 0; i < strips.strip_count; i++) {
        total_triangles += strips.strips[i].index_count - 2;
        total_indices += strips.strips[i].index_count;
    }

    TEST_ASSERT(total_triangles == mesh.triangle_count, "All triangles should be in strips");

    /* Calculate efficiency */
    size_t original_indices = mesh.triangle_count * 3;
    printf("  Original indices: %zu\n", original_indices);
    printf("  Strip indices: %zu\n", total_indices);
    printf("  Savings: %.2f%%\n",
           100.0f * (original_indices - total_indices) / original_indices);

    TEST_ASSERT(total_indices <= original_indices,
                "Strips should not use more indices than triangles");

    geomlib_tristrip_group_free(&strips);
    geomlib_trimesh_free(&mesh);

    printf("  PASS\n");
    return 0;
}

static int test_reorder_vertices() {
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 12, 2);
    for (int i = 0; i < 6; i++) {
        geomlib_vec3 normal = {float(i + 1), 0, 0};
        geomlib_vec2 uv = {float(i + 2), 0};
        geomlib_trimesh_add_vertex(&mesh, {float(i), 0, 0}, &normal, &uv);
    }
    geomlib_trimesh_add_triangle(&mesh, 2, 3, 4);
    geomlib_trimesh_add_triangle(&mesh, 4, 1, 2);
    const uint32_t bad[] = {1, 1}, order[] = {1, 0};
    TEST_ASSERT(geomlib_reorder_vertices(&mesh, bad) == GEOMLIB_MESH_ERROR_INVALID &&
                mesh.vertices[0].x == 0 && mesh.indices[0] == 2,
                "Invalid triangle permutation should leave the mesh unchanged");
    TEST_ASSERT(geomlib_reorder_vertices(&mesh, order) == GEOMLIB_MESH_SUCCESS,
                "Reorder vertices by the requested traversal");
    const int expected[] = {4, 1, 2, 3, 0, 5};
    for (int i = 0; i < 6; i++) {
        TEST_ASSERT(mesh.vertices[i].x == expected[i] && mesh.normals[i].x == expected[i] + 1 &&
                    mesh.texcoords[i].x == expected[i] + 2, "Vertex attributes remain aligned");
    }
    const uint32_t corners[] = {2, 3, 4, 4, 1, 2};
    for (size_t i = 0; i < 6; i++)
        TEST_ASSERT(mesh.vertices[mesh.indices[i]].x == corners[i], "Original triangle order and geometry remain intact");
    TEST_ASSERT(geomlib_optimize_vertex_fetch(&mesh) == GEOMLIB_MESH_SUCCESS,
                "Default vertex fetch uses the same implementation");
    TEST_ASSERT(geomlib_trimesh_add_vertex(&mesh, {}, nullptr, nullptr) == GEOMLIB_MESH_SUCCESS,
                "Reserved vertex storage remains usable after reordering");
    geomlib_trimesh_free(&mesh);
    return 0;
}

int main() {
    printf("=== Mesh Optimization Tests ===\n\n");

    int failures = 0;

    failures += test_triangle_strips();
    failures += test_index_compression();
    failures += test_reorder_vertices();
    failures += test_compression_sphere();
    failures += test_strips_sphere();

    printf("\n");
    if (failures == 0) {
        printf("=== ALL TESTS PASSED ===\n");
    } else {
        printf("=== %d TESTS FAILED ===\n", failures);
    }

    return failures;
}
