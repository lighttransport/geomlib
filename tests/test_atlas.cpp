/* Atlas statistics: normalized stretch is one for isometric UVs and sqrt(2) for a 2x stretch. */
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "geomlib/geomlib.h"

#define CHECK(condition, message) do { \
    if (!(condition)) { std::printf("FAIL: %s (line %d)\n", message, __LINE__); return 1; } \
} while (0)

static geomlib_trimesh make_square(float u_scale, float v_scale) {
    geomlib_trimesh mesh = {};
    geomlib_trimesh_init(&mesh, 4, 2);
    const geomlib_vec3 positions[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    for (int i = 0; i < 4; i++) {
        geomlib_vec2 uv = {positions[i].x * u_scale, positions[i].y * v_scale};
        geomlib_trimesh_add_vertex(&mesh, positions[i], nullptr, &uv);
    }
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 0, 2, 3);
    return mesh;
}

static geomlib_atlas make_atlas_for_all_triangles(const geomlib_trimesh *mesh) {
    geomlib_atlas atlas = {};
    atlas.chart_count = 1;
    atlas.charts = (geomlib_atlas_chart*)calloc(1, sizeof(geomlib_atlas_chart));
    atlas.charts[0].triangles = (uint32_t*)malloc(mesh->triangle_count * sizeof(uint32_t));
    atlas.charts[0].triangle_count = mesh->triangle_count;
    for (size_t i = 0; i < mesh->triangle_count; i++) atlas.charts[0].triangles[i] = (uint32_t)i;
    atlas.charts[0].uv_min = {0, 0};
    atlas.charts[0].uv_max = {1, 1};
    return atlas;
}

int main() {
    std::printf("=== Atlas Test Suite ===\n\n");

    geomlib_trimesh isometric = make_square(1.0f, 1.0f);
    geomlib_atlas atlas = make_atlas_for_all_triangles(&isometric);
    geomlib_atlas_compute_stats(&isometric, &atlas);
    CHECK(std::fabs(atlas.avg_stretch - 1.0f) < 1e-4f &&
          std::fabs(atlas.max_stretch - 1.0f) < 1e-4f,
          "Isometric UVs have unit stretch");
    geomlib_atlas_free(&atlas);
    geomlib_trimesh_free(&isometric);

    geomlib_trimesh stretched = make_square(2.0f, 1.0f);
    atlas = make_atlas_for_all_triangles(&stretched);
    geomlib_atlas_compute_stats(&stretched, &atlas);
    CHECK(std::fabs(atlas.max_stretch - std::sqrt(2.0f)) < 1e-3f,
          "A 2x UV stretch reports sqrt(2) maximum stretch");
    CHECK(atlas.avg_stretch > 1.05f && atlas.avg_stretch < 1.12f,
          "Average stretch reflects the anisotropic mapping");
    std::printf("  avg stretch %.4f, max %.4f\n", atlas.avg_stretch, atlas.max_stretch);
    geomlib_atlas_free(&atlas);
    geomlib_trimesh_free(&stretched);

    std::printf("\n=== ALL TESTS PASSED ===\n");
    return 0;
}
