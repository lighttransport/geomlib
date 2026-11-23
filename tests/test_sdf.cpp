/*
 * SDF generation modes: fast sweeping, sparse narrow band, feature extraction.
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

static double box_signed_distance(double x, double y, double z) {
    double q[3] = {fabs(x) - 0.5, fabs(y) - 0.5, fabs(z) - 0.5};
    double outside = 0.0;
    for (int i = 0; i < 3; i++) if (q[i] > 0.0) outside += q[i] * q[i];
    outside = sqrt(outside);
    double inside = fmin(fmax(q[0], fmax(q[1], q[2])), 0.0);
    return outside + inside;
}

static double mesh_volume(const geomlib_trimesh *mesh) {
    double volume = 0.0;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 a, b, c;
        geomlib_trimesh_get_triangle_vertices(mesh, i, &a, &b, &c);
        volume += (double)geomlib_vec3_dot(a, geomlib_vec3_cross(b, c)) / 6.0;
    }
    return fabs(volume);
}

static int test_sparse_band() {
    printf("Testing sparse narrow-band storage...\n");
    geomlib_trimesh sphere;
    TEST_ASSERT(geomlib_create_sphere(&sphere, 1.0f, 24, 12) == GEOMLIB_MESH_SUCCESS,
                "Sphere fixture should build");

    geomlib_sdf_params params = geomlib_sdf_params_default();
    params.resolution = 32;
    params.voxel_size = 0.0f;
    params.narrow_band = true;
    params.band_width = 3.0f;
    params.sparse = true;
    geomlib_sdf sparse = geomlib_mesh_to_sdf(&sphere, &params);
    TEST_ASSERT(sparse.is_sparse && sparse.values == nullptr && sparse.band != nullptr,
                "Sparse generation should omit the dense array");
    TEST_ASSERT(sparse.stored_voxels > 0 && sparse.stored_voxels < sparse.total_voxels,
                "Sparse storage should keep only the band");
    TEST_ASSERT(geomlib_sdf_stored_voxels(&sparse) == sparse.stored_voxels,
                "Stored-voxel accessor should match the grid");
    TEST_ASSERT(geomlib_sdf_is_defined(&sparse, 0, 0, 0) == false,
                "The padded corner is outside the band");

    geomlib_trimesh extracted;
    TEST_ASSERT(geomlib_sdf_extract_surface(&sparse, &extracted) == GEOMLIB_MESH_SUCCESS &&
                extracted.triangle_count > 0,
                "Sparse extraction should produce a surface");
    geomlib_manifold_analysis analysis = geomlib_analyze_manifold(&extracted);
    TEST_ASSERT(analysis.is_closed && analysis.is_manifold && analysis.is_oriented,
                "Sparse extraction should be watertight and oriented");
    geomlib_manifold_analysis_free(&analysis);

    geomlib_sdf_stats stats = geomlib_compute_sdf_stats(&sparse);
    TEST_ASSERT(stats.memory_bytes < sparse.total_voxels * sizeof(float),
                "Sparse stats should report less memory than dense");

    geomlib_trimesh_free(&extracted);
    geomlib_sdf_free(&sparse);
    geomlib_trimesh_free(&sphere);
    return 0;
}

static int test_fast_sweeping_accuracy() {
    printf("Testing fast sweeping against brute force...\n");
    geomlib_trimesh sphere;
    TEST_ASSERT(geomlib_create_sphere(&sphere, 1.0f, 24, 12) == GEOMLIB_MESH_SUCCESS,
                "Sphere fixture should build");

    geomlib_sdf_params params = geomlib_sdf_params_default();
    params.resolution = 24;
    params.voxel_size = 0.0f;
    geomlib_sdf brute = geomlib_mesh_to_sdf(&sphere, &params);
    params.algorithm = GEOMLIB_SDF_ALGORITHM_FAST_SWEEPING;
    geomlib_sdf sweep = geomlib_mesh_to_sdf(&sphere, &params);
    TEST_ASSERT(brute.values && sweep.values && brute.total_voxels == sweep.total_voxels,
                "Both modes should produce dense grids");

    float h = brute.voxel_size;
    double total_error = 0.0;
    double max_error = 0.0;
    size_t sign_matches = 0;
    for (size_t i = 0; i < brute.total_voxels; i++) {
        double error = fabs((double)brute.values[i] - (double)sweep.values[i]);
        total_error += error;
        if (error > max_error) max_error = error;
        if ((brute.values[i] < 0.0f) == (sweep.values[i] < 0.0f)) sign_matches++;
    }
    double mean_error = total_error / (double)brute.total_voxels;
    TEST_ASSERT(mean_error < 0.5 * h, "Mean sweeping error should be sub-voxel");
    TEST_ASSERT(max_error < 3.0 * h, "Max sweeping error should stay within a few voxels");
    TEST_ASSERT((double)sign_matches / (double)brute.total_voxels > 0.995,
                "Sweeping signs should agree with brute force almost everywhere");

    geomlib_sdf_free(&brute);
    geomlib_sdf_free(&sweep);
    geomlib_trimesh_free(&sphere);
    return 0;
}

static int test_sign_acceleration() {
    printf("Testing BVH-accelerated inside test...\n");
    geomlib_trimesh sphere;
    TEST_ASSERT(geomlib_create_sphere(&sphere, 1.0f, 16, 8) == GEOMLIB_MESH_SUCCESS,
                "Sphere fixture should build");
    geomlib_bvh bvh;
    TEST_ASSERT(geomlib_bvh_build(&bvh, &sphere, GEOMLIB_BVH_QUALITY_MEDIUM) ==
                GEOMLIB_MESH_SUCCESS, "BVH should build");

    const geomlib_vec3 points[] = {
        {0, 0, 0}, {0.5f, 0.2f, 0.1f}, {0.8f, 0.8f, 0.0f}, {2.0f, 0.0f, 0.0f},
        {-0.3f, -0.9f, 0.1f}, {0.0f, 0.0f, 1.5f}, {0.1f, -0.2f, 0.95f}
    };
    for (const geomlib_vec3 &point : points) {
        bool reference = geomlib_is_point_inside_mesh(&sphere, nullptr, point);
        bool accelerated = geomlib_is_point_inside_mesh(&sphere, &bvh, point);
        TEST_ASSERT(reference == accelerated, "BVH and brute-force parity must agree");
    }
    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&sphere);
    return 0;
}

static int test_feature_extraction() {
    printf("Testing sharp-feature-preserving extraction...\n");
    geomlib_trimesh box;
    TEST_ASSERT(geomlib_create_box(&box, 1.0f, 1.0f, 1.0f) == GEOMLIB_MESH_SUCCESS,
                "Box fixture should build");

    geomlib_sdf_params params = geomlib_sdf_params_default();
    params.resolution = 24;
    params.voxel_size = 0.0f;
    geomlib_sdf sdf = geomlib_mesh_to_sdf(&box, &params);
    TEST_ASSERT(sdf.values != nullptr, "Box SDF should generate");

    geomlib_trimesh baseline;
    TEST_ASSERT(geomlib_sdf_extract_surface(&sdf, &baseline) == GEOMLIB_MESH_SUCCESS,
                "Baseline extraction should succeed");

    geomlib_trimesh featured;
    TEST_ASSERT(geomlib_sdf_extract_surface_features(&sdf, &box, 45.0f, true, &featured) ==
                GEOMLIB_MESH_SUCCESS && featured.triangle_count > 0,
                "Feature extraction should succeed");

    double max_deviation = 0.0;
    for (size_t i = 0; i < featured.vertex_count; i++) {
        const geomlib_vec3 &v = featured.vertices[i];
        double deviation = fabs(box_signed_distance(v.x, v.y, v.z));
        if (deviation > max_deviation) max_deviation = deviation;
    }
    TEST_ASSERT(max_deviation < 2e-3,
                "Feature vertices should lie on the reference box surface");

    int corners_found = 0;
    for (int c = 0; c < 8; c++) {
        geomlib_vec3 corner = {(c & 1) ? 0.5f : -0.5f, (c & 2) ? 0.5f : -0.5f,
                               (c & 4) ? 0.5f : -0.5f};
        double best = 1e9;
        for (size_t i = 0; i < featured.vertex_count; i++) {
            double d = geomlib_vec3_length(geomlib_vec3_sub(featured.vertices[i], corner));
            if (d < best) best = d;
        }
        if (best < 0.03) corners_found++;
    }
    TEST_ASSERT(corners_found == 8, "Every sharp box corner should be surrounded by vertices");

    int edge_vertices = 0;
    for (size_t i = 0; i < featured.vertex_count; i++) {
        const geomlib_vec3 &v = featured.vertices[i];
        int coordinates = (fabs(fabs(v.x) - 0.5f) < 2e-3f) + (fabs(fabs(v.y) - 0.5f) < 2e-3f) +
                          (fabs(fabs(v.z) - 0.5f) < 2e-3f);
        if (coordinates >= 2) edge_vertices++;
    }
    TEST_ASSERT(edge_vertices >= 12, "Extracted vertices should snap onto box sharp edges");

    double volume = mesh_volume(&featured);
    TEST_ASSERT(volume > 0.97 && volume < 1.03, "Feature mesh should preserve box volume");

    geomlib_manifold_analysis analysis = geomlib_analyze_manifold(&featured);
    TEST_ASSERT(analysis.is_closed && analysis.is_manifold && analysis.is_oriented,
                "Feature mesh should stay watertight");
    geomlib_manifold_analysis_free(&analysis);

    geomlib_trimesh_free(&baseline);
    geomlib_trimesh_free(&featured);
    geomlib_sdf_free(&sdf);
    geomlib_trimesh_free(&box);
    return 0;
}

static int test_sparse_round_trip() {
    printf("Testing sparse SDF2 round trip...\n");
    geomlib_trimesh sphere;
    TEST_ASSERT(geomlib_create_sphere(&sphere, 1.0f, 16, 8) == GEOMLIB_MESH_SUCCESS,
                "Sphere fixture should build");
    geomlib_sdf_params params = geomlib_sdf_params_default();
    params.resolution = 20;
    params.voxel_size = 0.0f;
    params.narrow_band = true;
    params.sparse = true;
    geomlib_sdf sparse = geomlib_mesh_to_sdf(&sphere, &params);
    TEST_ASSERT(sparse.is_sparse, "Fixture should be sparse");

    TEST_ASSERT(geomlib_sdf_save(&sparse, "/tmp/geomlib_sparse_sdf.bin") ==
                GEOMLIB_MESH_SUCCESS, "Sparse save should succeed");
    geomlib_sdf loaded = {};
    TEST_ASSERT(geomlib_sdf_load(&loaded, "/tmp/geomlib_sparse_sdf.bin") ==
                GEOMLIB_MESH_SUCCESS, "Sparse load should succeed");
    TEST_ASSERT(loaded.is_sparse && loaded.values == nullptr &&
                loaded.stored_voxels == sparse.stored_voxels,
                "Loaded sparse grid should preserve band storage");
    for (uint32_t z = 0; z < sparse.resolution[2]; z++) {
        for (uint32_t y = 0; y < sparse.resolution[1]; y++) {
            for (uint32_t x = 0; x < sparse.resolution[0]; x++) {
                float a, b;
                bool da = geomlib_sdf_is_defined(&sparse, x, y, z);
                bool db = geomlib_sdf_is_defined(&loaded, x, y, z);
                TEST_ASSERT(da == db, "Defined set should survive the round trip");
                if (da) {
                    a = geomlib_sdf_query_distance_voxel(&sparse, x, y, z);
                    b = geomlib_sdf_query_distance_voxel(&loaded, x, y, z);
                    TEST_ASSERT(a == b, "Sample values should survive the round trip");
                }
            }
        }
    }

    TEST_ASSERT(geomlib_sdf_make_dense(&loaded) == GEOMLIB_MESH_SUCCESS &&
                !loaded.is_sparse && loaded.values != nullptr,
                "Materializing a band should produce a dense grid");

    geomlib_sdf_free(&loaded);
    geomlib_sdf_free(&sparse);
    geomlib_trimesh_free(&sphere);
    return 0;
}

int main() {
    printf("=== SDF Mode Test Suite ===\n\n");
    int failures = 0;
    failures += test_sparse_band();
    failures += test_fast_sweeping_accuracy();
    failures += test_sign_acceleration();
    failures += test_feature_extraction();
    failures += test_sparse_round_trip();
    printf("\n");
    if (!failures) printf("=== ALL TESTS PASSED ===\n");
    else printf("=== %d TESTS FAILED ===\n", failures);
    return failures;
}
