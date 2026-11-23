/* Smoothing: Taubin volume preservation, HC correction, feature/boundary locks. */
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include "geomlib/geomlib.h"

#define CHECK(condition, message) do { \
    if (!(condition)) { std::printf("FAIL: %s (line %d)\n", message, __LINE__); return 1; } \
} while (0)

static double mesh_volume(const geomlib_trimesh *mesh) {
    double volume = 0.0;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 a, b, c;
        geomlib_trimesh_get_triangle_vertices(mesh, i, &a, &b, &c);
        volume += geomlib_vec3_dot(a, geomlib_vec3_cross(b, c)) / 6.0;
    }
    return std::fabs(volume);
}

static void weld_mesh(geomlib_trimesh *mesh) {
    geomlib_weld_result result = geomlib_weld_vertices(mesh, 1e-5f);
    geomlib_apply_weld_result(mesh, &result);
    geomlib_weld_result_free(&result);
}

static geomlib_trimesh make_sphere() {
    geomlib_trimesh mesh = {};
    geomlib_create_sphere(&mesh, 1.0f, 16, 12);
    weld_mesh(&mesh);
    return mesh;
}

static int test_taubin_volume() {
    geomlib_trimesh mesh = make_sphere();
    double before = mesh_volume(&mesh);
    auto params = geomlib_smoothing_params_default(GEOMLIB_SMOOTH_TAUBIN);
    params.iterations = 20;
    CHECK(geomlib_smooth_mesh(&mesh, &params) == GEOMLIB_MESH_SUCCESS, "Taubin should run");
    double after = mesh_volume(&mesh);
    CHECK(std::fabs(after - before) / before < 0.03, "Taubin preserves volume on a sphere");
    geomlib_trimesh_free(&mesh);

    mesh = make_sphere();
    before = mesh_volume(&mesh);
    params = geomlib_smoothing_params_default(GEOMLIB_SMOOTH_LAPLACIAN);
    params.iterations = 20;
    CHECK(geomlib_smooth_mesh(&mesh, &params) == GEOMLIB_MESH_SUCCESS, "Laplacian should run");
    after = mesh_volume(&mesh);
    CHECK(after < before * 0.97, "Plain Laplacian shrinks the sphere");
    geomlib_trimesh_free(&mesh);
    return 0;
}

/* A welded box has eight vertices, all on sharp edges. */
static int test_feature_lock() {
    geomlib_trimesh mesh = {};
    geomlib_create_box(&mesh, 1.0f, 1.0f, 1.0f);
    weld_mesh(&mesh);
    CHECK(mesh.vertex_count == 8, "The box welds to eight corners");
    std::vector<geomlib_vec3> original(mesh.vertices, mesh.vertices + mesh.vertex_count);

    auto params = geomlib_smoothing_params_default(GEOMLIB_SMOOTH_LAPLACIAN);
    params.iterations = 5;
    params.preserve_boundaries = false;
    params.preserve_features = true;
    params.feature_angle = 30.0f;
    CHECK(geomlib_smooth_mesh(&mesh, &params) == GEOMLIB_MESH_SUCCESS, "Feature-locked smoothing runs");
    for (size_t i = 0; i < mesh.vertex_count; i++)
        CHECK(std::memcmp(&original[i], &mesh.vertices[i], sizeof(geomlib_vec3)) == 0,
              "Sharp-edge vertices stay fixed");
    geomlib_trimesh_free(&mesh);

    mesh = {};
    geomlib_create_box(&mesh, 1.0f, 1.0f, 1.0f);
    weld_mesh(&mesh);
    params.preserve_features = false;
    CHECK(geomlib_smooth_mesh(&mesh, &params) == GEOMLIB_MESH_SUCCESS, "Free smoothing runs");
    bool moved = false;
    for (size_t i = 0; i < mesh.vertex_count; i++)
        moved |= std::fabs(mesh.vertices[i].x) != 0.5f ||
                 std::fabs(mesh.vertices[i].y) != 0.5f ||
                 std::fabs(mesh.vertices[i].z) != 0.5f;
    CHECK(moved, "Vertices move when feature preservation is off");
    geomlib_trimesh_free(&mesh);
    return 0;
}

static geomlib_trimesh make_plane(int n) {
    geomlib_trimesh mesh = {};
    geomlib_trimesh_init(&mesh, n * n, (n - 1) * (n - 1) * 2);
    for (int y = 0; y < n; y++)
        for (int x = 0; x < n; x++)
            geomlib_trimesh_add_vertex(&mesh, {(float)x, (float)y, 0.05f * x * y}, nullptr, nullptr);
    for (int y = 0; y + 1 < n; y++)
        for (int x = 0; x + 1 < n; x++) {
            uint32_t a = y * n + x, b = a + 1, c = a + n, d = c + 1;
            geomlib_trimesh_add_triangle(&mesh, a, b, d);
            geomlib_trimesh_add_triangle(&mesh, a, d, c);
        }
    return mesh;
}

static int test_boundary_lock() {
    geomlib_trimesh mesh = make_plane(5);
    std::vector<geomlib_vec3> original(mesh.vertices, mesh.vertices + mesh.vertex_count);
    auto params = geomlib_smoothing_params_default(GEOMLIB_SMOOTH_LAPLACIAN);
    params.iterations = 3;
    params.preserve_boundaries = true;
    CHECK(geomlib_smooth_mesh(&mesh, &params) == GEOMLIB_MESH_SUCCESS, "Plane smoothing runs");
    for (int y = 0; y < 5; y++)
        for (int x = 0; x < 5; x++) {
            size_t i = y * 5 + x;
            if (x == 0 || y == 0 || x == 4 || y == 4)
                CHECK(std::memcmp(&original[i], &mesh.vertices[i], sizeof(geomlib_vec3)) == 0,
                      "Boundary vertices stay fixed");
            else
                CHECK(std::memcmp(&original[i], &mesh.vertices[i], sizeof(geomlib_vec3)) != 0,
                      "Interior vertices move");
        }
    CHECK(geomlib_smooth_mesh(&mesh, nullptr) == GEOMLIB_MESH_ERROR_INVALID,
          "Null parameters are rejected");
    params.iterations = -1;
    CHECK(geomlib_smooth_mesh(&mesh, &params) == GEOMLIB_MESH_ERROR_INVALID,
          "Negative iterations are rejected");
    geomlib_trimesh_free(&mesh);
    return 0;
}

int main() {
    std::printf("=== Smoothing Test Suite ===\n\n");
    int failures = 0;
    failures += test_taubin_volume();
    failures += test_feature_lock();
    failures += test_boundary_lock();
    if (!failures) std::printf("=== ALL TESTS PASSED ===\n");
    else std::printf("=== %d TESTS FAILED ===\n", failures);
    return failures;
}
