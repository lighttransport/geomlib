/* QEM decimation: edge collapses, boundary preservation, error limits. */
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include "geomlib/geomlib.h"

#define CHECK(condition, message) do { \
    if (!(condition)) { std::printf("FAIL: %s (line %d)\n", message, __LINE__); return 1; } \
} while (0)

static void weld_mesh(geomlib_trimesh *mesh) {
    geomlib_weld_result result = geomlib_weld_vertices(mesh, 1e-5f);
    geomlib_apply_weld_result(mesh, &result);
    geomlib_weld_result_free(&result);
}

static double mesh_volume(const geomlib_trimesh *mesh) {
    double volume = 0.0;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 a, b, c;
        geomlib_trimesh_get_triangle_vertices(mesh, i, &a, &b, &c);
        volume += geomlib_vec3_dot(a, geomlib_vec3_cross(b, c)) / 6.0;
    }
    return std::fabs(volume);
}

static bool edge_manifold(const geomlib_trimesh *mesh) {
    std::map<uint64_t, int> edges;
    std::set<std::array<uint32_t, 3>> faces;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t v[3] = {mesh->indices[i * 3], mesh->indices[i * 3 + 1], mesh->indices[i * 3 + 2]};
        std::array<uint32_t, 3> key = {v[0], v[1], v[2]};
        std::sort(key.begin(), key.end());
        if (!faces.insert(key).second) return false;
        for (int j = 0; j < 3; j++) {
            uint32_t a = v[j], b = v[(j + 1) % 3];
            if (a > b) std::swap(a, b);
            edges[((uint64_t)a << 32) | b]++;
        }
    }
    for (const auto &entry : edges) if (entry.second > 2) return false;
    return true;
}

static geomlib_trimesh make_closed_sphere() {
    geomlib_trimesh mesh = {};
    geomlib_trimesh_init(&mesh, 6, 8);
    for (geomlib_vec3 p : {geomlib_vec3{1, 0, 0}, geomlib_vec3{0, 1, 0},
                           geomlib_vec3{-1, 0, 0}, geomlib_vec3{0, -1, 0},
                           geomlib_vec3{0, 0, 1}, geomlib_vec3{0, 0, -1}})
        geomlib_trimesh_add_vertex(&mesh, p, nullptr, nullptr);
    for (uint32_t i = 0; i < 4; i++) {
        geomlib_trimesh_add_triangle(&mesh, i, (i + 1) % 4, 4);
        geomlib_trimesh_add_triangle(&mesh, (i + 1) % 4, i, 5);
    }
    for (int level = 0; level < 4; level++) {
        geomlib_trimesh next = {};
        if (geomlib_subdivide_loop(&mesh, &next) != GEOMLIB_MESH_SUCCESS) {
            geomlib_trimesh_free(&mesh);
            return {};
        }
        geomlib_trimesh_free(&mesh);
        mesh = next;
    }
    for (size_t i = 0; i < mesh.vertex_count; i++) {
        float length = geomlib_vec3_length(mesh.vertices[i]);
        if (length > 1e-6f) mesh.vertices[i] = geomlib_vec3_mul(mesh.vertices[i], 1.0f / length);
    }
    return mesh;
}

static int test_sphere_decimation() {
    geomlib_trimesh mesh = make_closed_sphere();
    size_t original = mesh.triangle_count;
    double original_volume = mesh_volume(&mesh);
    CHECK(original > 900, "Sphere fixture has enough triangles");

    geomlib_decimation_stats stats;
    geomlib_decimation_params params = geomlib_decimation_params_default(0.0f);
    params.target_triangle_count = 120;
    CHECK(geomlib_decimate_qem(&mesh, &params, &stats) == GEOMLIB_MESH_SUCCESS,
          "QEM decimation succeeds");
    CHECK(mesh.triangle_count <= 130 && mesh.triangle_count > 0,
          "Decimation reaches the target triangle count");
    CHECK(stats.edges_collapsed > 0 && stats.final_triangles == mesh.triangle_count &&
          stats.original_triangles == original, "Statistics track the reduction");
    CHECK(edge_manifold(&mesh), "Decimated mesh stays free of duplicate faces");
    geomlib_manifold_analysis analysis = geomlib_analyze_manifold(&mesh);
    CHECK(analysis.is_closed && analysis.is_manifold, "Closed sphere stays closed and manifold");
    geomlib_manifold_analysis_free(&analysis);
    CHECK(std::fabs(mesh_volume(&mesh) - original_volume) < 0.2 * original_volume,
          "Decimation keeps the sphere volume within 20 percent");
    geomlib_trimesh_free(&mesh);
    return 0;
}

static int test_boundary_preservation() {
    geomlib_trimesh mesh = {};
    geomlib_create_plane(&mesh, 4.0f, 4.0f, 6);
    geomlib_boundary_loops before = geomlib_extract_boundary_loops(&mesh);
    size_t before_vertices = before.loop_count ? before.loops[0].vertex_count : 0;
    geomlib_boundary_loops_free(&before);
    CHECK(before_vertices >= 4, "Plane has a boundary loop");

    geomlib_decimation_params params = geomlib_decimation_params_default(0.0f);
    params.target_triangle_count = mesh.triangle_count / 2;
    params.preserve_boundaries = true;
    CHECK(geomlib_decimate_qem(&mesh, &params, NULL) == GEOMLIB_MESH_SUCCESS,
          "Boundary-preserving decimation succeeds");
    geomlib_boundary_loops after = geomlib_extract_boundary_loops(&mesh);
    CHECK(after.loop_count == 1 && after.loops[0].vertex_count == before_vertices,
          "Boundary loop is unchanged when preserving boundaries");
    /* The four corners must still exist. */
    for (int c = 0; c < 4; c++) {
        geomlib_vec3 corner = {(c & 1) ? 2.0f : -2.0f, 0.0f, (c & 2) ? 2.0f : -2.0f};
        bool found = false;
        for (size_t i = 0; i < mesh.vertex_count; i++)
            found |= geomlib_vec3_length(geomlib_vec3_sub(mesh.vertices[i], corner)) < 1e-4f;
        CHECK(found, "Boundary corners survive preserved decimation");
    }
    geomlib_boundary_loops_free(&after);
    geomlib_trimesh_free(&mesh);
    return 0;
}

static int test_error_limit() {
    geomlib_trimesh mesh = {};
    geomlib_create_sphere(&mesh, 1.0f, 16, 8);
    size_t original = mesh.triangle_count;
    geomlib_decimation_params params = geomlib_decimation_params_default(0.0f);
    params.target_triangle_count = original / 4;
    params.max_error = 1e-8f;
    CHECK(geomlib_decimate_qem(&mesh, &params, NULL) == GEOMLIB_MESH_SUCCESS,
          "Error-limited decimation succeeds");
    CHECK(mesh.triangle_count >= original / 2,
          "A very small error limit stops before aggressive reduction");
    geomlib_trimesh_free(&mesh);
    return 0;
}

int main() {
    std::printf("=== Decimation Test Suite ===\n\n");
    int failures = 0;
    failures += test_sphere_decimation();
    failures += test_boundary_preservation();
    failures += test_error_limit();
    if (!failures) std::printf("=== ALL TESTS PASSED ===\n");
    else std::printf("=== %d TESTS FAILED ===\n", failures);
    return failures;
}
