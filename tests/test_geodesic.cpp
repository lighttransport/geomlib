/* Geodesics: heat-method distances against analytic plane and sphere metrics. */
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include "geomlib/geomlib.h"

#define CHECK(condition, message) do { \
    if (!(condition)) { std::printf("FAIL: %s (line %d)\n", message, __LINE__); return 1; } \
} while (0)

static geomlib_trimesh make_plane(int n) {
    geomlib_trimesh mesh = {};
    geomlib_trimesh_init(&mesh, n * n, (n - 1) * (n - 1) * 2);
    for (int y = 0; y < n; y++)
        for (int x = 0; x < n; x++)
            geomlib_trimesh_add_vertex(&mesh, {(float)x, (float)y, 0.0f}, nullptr, nullptr);
    for (int y = 0; y + 1 < n; y++)
        for (int x = 0; x + 1 < n; x++) {
            uint32_t a = y * n + x, b = a + 1, c = a + n, d = c + 1;
            geomlib_trimesh_add_triangle(&mesh, a, b, d);
            geomlib_trimesh_add_triangle(&mesh, a, d, c);
        }
    return mesh;
}

static int test_plane() {
    const int n = 11;
    geomlib_trimesh mesh = make_plane(n);
    uint32_t source = (uint32_t)(5 * n + 5);
    geomlib_geodesic_distances geo = geomlib_geodesic_create(mesh.vertex_count);
    CHECK(geomlib_compute_geodesic_from_vertex(&mesh, source, &geo) == GEOMLIB_MESH_SUCCESS,
          "Heat geodesics should succeed");
    CHECK(geo.source_count == 1 && geo.source_vertices[0] == source,
          "Source metadata is stored");
    CHECK(geo.distances[source] < 1e-4f, "The source distance is zero");
    double total_error = 0.0;
    int samples = 0;
    for (size_t i = 0; i < mesh.vertex_count; i++) {
        double dx = mesh.vertices[i].x - mesh.vertices[source].x;
        double dy = mesh.vertices[i].y - mesh.vertices[source].y;
        double exact = std::sqrt(dx * dx + dy * dy);
        if (exact < 1.0) continue;
        double error = std::fabs(geo.distances[i] - exact) / exact;
        total_error += error;
        samples++;
    }
    CHECK(samples > 50 && total_error / samples < 0.08,
          "Mean plane error stays below 8 percent");
    geomlib_geodesic_free(&geo);
    geomlib_trimesh_free(&mesh);
    return 0;
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
    for (int level = 0; level < 3; level++) {
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

static int test_sphere() {
    geomlib_trimesh mesh = make_closed_sphere();
    CHECK(mesh.vertex_count > 100 && mesh.triangle_count > 0, "Sphere fixture builds");
    geomlib_manifold_analysis analysis = geomlib_analyze_manifold(&mesh);
    CHECK(analysis.is_closed, "The sphere fixture is closed");
    geomlib_manifold_analysis_free(&analysis);

    uint32_t source = 0;
    for (size_t i = 1; i < mesh.vertex_count; i++)
        if (mesh.vertices[i].z > mesh.vertices[source].z) source = (uint32_t)i;
    geomlib_geodesic_distances geo = geomlib_geodesic_create(mesh.vertex_count);
    CHECK(geomlib_compute_geodesic_from_vertex(&mesh, source, &geo) == GEOMLIB_MESH_SUCCESS,
          "Sphere geodesics should succeed");
    double total_error = 0.0, max_error = 0.0;
    int samples = 0;
    for (size_t i = 0; i < mesh.vertex_count; i++) {
        double dot = geomlib_vec3_dot(mesh.vertices[source], mesh.vertices[i]);
        dot = std::fmax(-1.0, std::fmin(1.0, dot));
        double exact = std::acos(dot);
        if (exact < 0.5) continue;
        double error = std::fabs(geo.distances[i] - exact) / exact;
        total_error += error;
        max_error = std::fmax(max_error, error);
        samples++;
    }
    CHECK(samples > 50 && total_error / samples < 0.08,
          "Mean sphere error stays below 8 percent");
    CHECK(max_error < 0.30, "Worst sphere error stays below 30 percent");
    geomlib_geodesic_free(&geo);
    geomlib_trimesh_free(&mesh);
    return 0;
}

static int test_sources_and_time_step() {
    geomlib_trimesh mesh = make_plane(7);
    uint32_t sources[2] = {8, 40};
    geomlib_geodesic_distances geo = geomlib_geodesic_create(mesh.vertex_count);
    CHECK(geomlib_compute_geodesic_heat(&mesh, sources, 2, 0.0f, &geo) == GEOMLIB_MESH_SUCCESS,
          "Multi-source heat geodesics should succeed");
    CHECK(geo.distances[sources[0]] < 1e-3f && geo.distances[sources[1]] < 1e-3f,
          "Both sources sit at zero distance");

    float mean = geomlib_compute_mean_edge_length_for_geodesic(&mesh);
    geomlib_geodesic_distances explicit_geo = geomlib_geodesic_create(mesh.vertex_count);
    CHECK(geomlib_compute_geodesic_heat(&mesh, sources, 2, mean * mean, &explicit_geo) ==
          GEOMLIB_MESH_SUCCESS, "Explicit time step should succeed");
    for (size_t i = 0; i < mesh.vertex_count; i++)
        CHECK(std::fabs(geo.distances[i] - explicit_geo.distances[i]) < 1e-3f,
              "Automatic time step matches the explicit mean-edge-length value");

    geomlib_geodesic_free(&geo);
    geomlib_geodesic_free(&explicit_geo);
    geomlib_trimesh_free(&mesh);
    return 0;
}

int main() {
    std::printf("=== Geodesic Test Suite ===\n\n");
    int failures = 0;
    failures += test_plane();
    failures += test_sphere();
    failures += test_sources_and_time_step();
    if (!failures) std::printf("=== ALL TESTS PASSED ===\n");
    else std::printf("=== %d TESTS FAILED ===\n", failures);
    return failures;
}
