/* Overdraw: rasterized fragment counts and front-to-back reordering. */
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include "geomlib/geomlib.h"

#define CHECK(condition, message) do { \
    if (!(condition)) { std::printf("FAIL: %s (line %d)\n", message, __LINE__); return 1; } \
} while (0)

static geomlib_trimesh make_layer(float z, bool facing_camera) {
    geomlib_trimesh mesh = {};
    geomlib_trimesh_init(&mesh, 3, 1);
    geomlib_trimesh_add_vertex(&mesh, {0, 0, z}, nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, {1, 0, z}, nullptr, nullptr);
    geomlib_trimesh_add_vertex(&mesh, {0, 1, z}, nullptr, nullptr);
    /* Camera is at negative z: normals must point toward -z. */
    if (facing_camera)
        geomlib_trimesh_add_triangle(&mesh, 0, 2, 1);
    else
        geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    return mesh;
}

static geomlib_trimesh make_two_layers() {
    geomlib_trimesh mesh = {};
    geomlib_trimesh_init(&mesh, 6, 2);
    for (float z : {0.0f, 1.0f}) {
        uint32_t base = (uint32_t)mesh.vertex_count;
        geomlib_trimesh_add_vertex(&mesh, {0, 0, z}, nullptr, nullptr);
        geomlib_trimesh_add_vertex(&mesh, {1, 0, z}, nullptr, nullptr);
        geomlib_trimesh_add_vertex(&mesh, {0, 1, z}, nullptr, nullptr);
        geomlib_trimesh_add_triangle(&mesh, base + 0, base + 2, base + 1);
    }
    return mesh;
}

int main() {
    std::printf("=== Overdraw Test Suite ===\n\n");
    geomlib_vec3 eye = {0.5f, 0.5f, -5.0f};

    geomlib_trimesh single = make_layer(0.0f, true);
    geomlib_overdraw_stats stats = geomlib_analyze_overdraw(&single, eye, 64);
    CHECK(stats.visible_pixels > 0, "A single layer covers pixels");
    CHECK(std::fabs(stats.avg_overdraw - 1.0f) < 1e-4f &&
          std::fabs(stats.max_overdraw - 1.0f) < 1e-4f,
          "A single layer has unit overdraw");
    geomlib_trimesh_free(&single);

    geomlib_trimesh layers = make_two_layers();
    stats = geomlib_analyze_overdraw(&layers, eye, 64);
    CHECK(stats.visible_pixels > 0 &&
          std::fabs(stats.avg_overdraw - 2.0f) < 1e-4f &&
          std::fabs(stats.max_overdraw - 2.0f) < 1e-4f,
          "Two identical layers have overdraw two");
    CHECK(stats.total_shaded_pixels == stats.visible_pixels * 2,
          "Fragment count matches two layers");

    CHECK(geomlib_optimize_overdraw(&layers, eye) == GEOMLIB_MESH_SUCCESS,
          "Front-to-back reordering succeeds");
    for (int j = 0; j < 3; j++) {
        CHECK(layers.vertices[layers.indices[j]].z == 0.0f,
              "Near layer is emitted first");
        CHECK(layers.vertices[layers.indices[3 + j]].z == 1.0f,
              "Far layer is emitted last");
    }
    geomlib_overdraw_stats reordered = geomlib_analyze_overdraw(&layers, eye, 64);
    CHECK(std::fabs(reordered.avg_overdraw - 2.0f) < 1e-4f,
          "Reordering does not change overdraw");
    geomlib_trimesh_free(&layers);

    geomlib_trimesh back = make_layer(0.0f, false);
    stats = geomlib_analyze_overdraw(&back, eye, 64);
    CHECK(stats.visible_pixels == 0 && stats.avg_overdraw == 0.0f,
          "Back-facing geometry is culled");
    geomlib_trimesh_free(&back);

    std::printf("\n=== ALL TESTS PASSED ===\n");
    return 0;
}
