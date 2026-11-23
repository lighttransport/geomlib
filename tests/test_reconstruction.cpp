/* Surface reconstruction: SDF marching cubes, alpha shape, empty-ball BPA, greedy projection. */
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <vector>
#include "geomlib/geomlib.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define CHECK(condition, message) do { \
    if (!(condition)) { std::printf("FAIL: %s (line %d)\n", message, __LINE__); return 1; } \
} while (0)

static uint64_t edge_key(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return ((uint64_t)a << 32) | b;
}

static bool manifold_faces(const geomlib_trimesh *mesh, bool require_closed) {
    std::map<uint64_t, int> edges;
    std::set<std::array<uint32_t, 3>> faces;
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        uint32_t v[3] = {mesh->indices[i * 3], mesh->indices[i * 3 + 1], mesh->indices[i * 3 + 2]};
        std::array<uint32_t, 3> key = {v[0], v[1], v[2]};
        std::sort(key.begin(), key.end());
        if (!faces.insert(key).second) return false;
        for (int j = 0; j < 3; j++) edges[edge_key(v[j], v[(j + 1) % 3])]++;
    }
    for (const auto &entry : edges) {
        if (entry.second > 2) return false;
        if (require_closed && entry.second != 2) return false;
    }
    return true;
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

static geomlib_pointcloud make_plane_cloud(int n) {
    geomlib_pointcloud pc;
    geomlib_pointcloud_init(&pc, n * n);
    geomlib_vec3 normal = {0, 0, 1};
    for (int y = 0; y < n; y++)
        for (int x = 0; x < n; x++)
            geomlib_pointcloud_add_point(&pc, {(float)x, (float)y, 0.0f}, &normal, 1.0f);
    return pc;
}

static geomlib_pointcloud make_sphere_cloud(size_t count) {
    geomlib_pointcloud pc;
    geomlib_pointcloud_init(&pc, count);
    for (size_t i = 0; i < count; i++) {
        double phi = std::acos(1.0 - 2.0 * ((double)i + 0.5) / (double)count);
        double theta = M_PI * (1.0 + std::sqrt(5.0)) * (double)i;
        double jitter = 1.0 + 0.02 * std::sin(13.0 * (double)i);
        geomlib_vec3 p = {(float)(jitter * std::sin(phi) * std::cos(theta)),
                          (float)(jitter * std::sin(phi) * std::sin(theta)),
                          (float)(jitter * std::cos(phi))};
        geomlib_vec3 normal = {(float)(std::sin(phi) * std::cos(theta)),
                               (float)(std::sin(phi) * std::sin(theta)),
                               (float)std::cos(phi)};
        geomlib_pointcloud_add_point(&pc, p, &normal, 1.0f);
    }
    return pc;
}

static bool faces_normal_to(const geomlib_trimesh *mesh, geomlib_vec3 expected) {
    for (size_t i = 0; i < mesh->triangle_count; i++) {
        geomlib_vec3 a, b, c;
        geomlib_trimesh_get_triangle_vertices(mesh, i, &a, &b, &c);
        geomlib_vec3 n = geomlib_vec3_cross(geomlib_vec3_sub(b, a), geomlib_vec3_sub(c, a));
        float length = geomlib_vec3_length(n);
        if (!(length > 1e-9f)) continue;
        if (std::fabs(geomlib_vec3_dot(geomlib_vec3_mul(n, 1.0f / length), expected)) < 0.99f)
            return false;
    }
    return true;
}

static int test_plane_methods() {
    geomlib_pointcloud pc = make_plane_cloud(6);
    geomlib_vec3 up = {0, 0, 1};

    geomlib_surface_params params = geomlib_surface_params_default(GEOMLIB_SURFACE_GREEDY_PROJECTION);
    geomlib_trimesh mesh = {};
    CHECK(geomlib_reconstruct_surface(&pc, &params, &mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.triangle_count > 0, "Greedy projection produces triangles");
    CHECK(manifold_faces(&mesh, false), "Greedy output is edge-manifold without duplicate faces");
    CHECK(faces_normal_to(&mesh, up), "Greedy triangles follow the plane normal");
    geomlib_trimesh_free(&mesh);

    params = geomlib_surface_params_default(GEOMLIB_SURFACE_BALL_PIVOTING);
    params.ball_radius = 1.5f;
    CHECK(geomlib_reconstruct_surface(&pc, &params, &mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.triangle_count > 0, "Ball pivoting produces triangles");
    CHECK(manifold_faces(&mesh, false), "Ball pivoting output is edge-manifold");
    CHECK(faces_normal_to(&mesh, up), "Ball pivoting triangles follow the plane normal");
    geomlib_trimesh_free(&mesh);

    params = geomlib_surface_params_default(GEOMLIB_SURFACE_MARCHING_CUBES);
    params.voxel_size = 0.25f;
    CHECK(geomlib_reconstruct_surface(&pc, &params, &mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.triangle_count > 0, "Marching cubes extracts a surface");
    CHECK(manifold_faces(&mesh, false), "Marching cubes output is manifold");
    geomlib_trimesh_free(&mesh);

    params = geomlib_surface_params_default(GEOMLIB_SURFACE_ALPHA_SHAPE);
    params.alpha = 2.0f;
    CHECK(geomlib_reconstruct_surface(&pc, &params, &mesh) == GEOMLIB_MESH_ERROR_DEGENERATE,
          "A coplanar cloud cannot form an alpha shape");
    geomlib_pointcloud_free(&pc);
    return 0;
}

static int test_sphere_methods() {
    geomlib_pointcloud pc = make_sphere_cloud(200);
    geomlib_trimesh mesh = {};

    geomlib_surface_params params = geomlib_surface_params_default(GEOMLIB_SURFACE_MARCHING_CUBES);
    CHECK(geomlib_reconstruct_surface(&pc, &params, &mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.triangle_count > 0, "Marching cubes reconstructs the sphere");
    CHECK(manifold_faces(&mesh, true), "Sphere marching cubes is closed and manifold");
    double volume = mesh_volume(&mesh);
    CHECK(volume > 1.5 && volume < 8.0, "Sphere marching-cubes volume is in a plausible range");
    geomlib_trimesh_free(&mesh);

    params = geomlib_surface_params_default(GEOMLIB_SURFACE_ALPHA_SHAPE);
    params.alpha = 1.0f;
    CHECK(geomlib_reconstruct_surface(&pc, &params, &mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.triangle_count > 0, "Alpha shape reconstructs the sphere");
    CHECK(manifold_faces(&mesh, false), "Alpha shape is edge-manifold");
    volume = mesh_volume(&mesh);
    CHECK(volume > 1.0 && volume < 8.0, "Alpha-shape volume is in a plausible range");
    geomlib_trimesh_free(&mesh);

    params = geomlib_surface_params_default(GEOMLIB_SURFACE_BALL_PIVOTING);
    params.ball_radius = 0.6f;
    CHECK(geomlib_reconstruct_surface(&pc, &params, &mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.triangle_count > 0, "Ball pivoting reconstructs the sphere");
    CHECK(manifold_faces(&mesh, false), "Sphere ball pivoting is edge-manifold");
    geomlib_trimesh_free(&mesh);

    params = geomlib_surface_params_default(GEOMLIB_SURFACE_GREEDY_PROJECTION);
    CHECK(geomlib_reconstruct_surface(&pc, &params, &mesh) == GEOMLIB_MESH_SUCCESS &&
          mesh.triangle_count > 0, "Greedy projection reconstructs the sphere");
    CHECK(manifold_faces(&mesh, false), "Sphere greedy projection is edge-manifold");
    geomlib_trimesh_free(&mesh);
    geomlib_pointcloud_free(&pc);
    return 0;
}

static int test_auto_target() {
    geomlib_pointcloud pc = make_sphere_cloud(200);
    geomlib_trimesh coarse = {}, fine = {};
    CHECK(geomlib_reconstruct_surface_auto(&pc, 100, &coarse) == GEOMLIB_MESH_SUCCESS &&
          coarse.triangle_count > 0, "Auto reconstruction with a small target succeeds");
    CHECK(geomlib_reconstruct_surface_auto(&pc, 800, &fine) == GEOMLIB_MESH_SUCCESS &&
          fine.triangle_count > 0, "Auto reconstruction with a large target succeeds");
    CHECK(coarse.triangle_count <= fine.triangle_count,
          "target_faces controls the reconstruction density");
    std::printf("  auto target 100 -> %zu triangles, 800 -> %zu triangles\n",
                coarse.triangle_count, fine.triangle_count);
    geomlib_trimesh_free(&coarse);
    geomlib_trimesh_free(&fine);
    geomlib_pointcloud_free(&pc);
    return 0;
}

static int test_invalid_inputs() {
    geomlib_pointcloud pc;
    geomlib_pointcloud_init(&pc, 0);
    geomlib_pointcloud_add_point(&pc, {0, 0, 0}, nullptr, 1.0f);
    geomlib_pointcloud_add_point(&pc, {1, 0, 0}, nullptr, 1.0f);
    geomlib_surface_params params =
        geomlib_surface_params_default(GEOMLIB_SURFACE_GREEDY_PROJECTION);
    geomlib_trimesh mesh = {};
    CHECK(geomlib_reconstruct_surface(&pc, &params, &mesh) == GEOMLIB_MESH_ERROR_INVALID,
          "Fewer than three points are rejected");
    geomlib_pointcloud_free(&pc);
    return 0;
}

int main() {
    std::printf("=== Reconstruction Test Suite ===\n\n");
    int failures = 0;
    failures += test_plane_methods();
    failures += test_sphere_methods();
    failures += test_auto_target();
    failures += test_invalid_inputs();
    if (!failures) std::printf("=== ALL TESTS PASSED ===\n");
    else std::printf("=== %d TESTS FAILED ===\n", failures);
    return failures;
}
