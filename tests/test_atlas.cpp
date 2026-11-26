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

    /* Wave 1: full atlas pipeline with all segmentation and packing modes. */
    {
        geomlib_trimesh bent;
        geomlib_trimesh_init(&bent, 6, 4);
        /* Two planar patches meeting at a right angle along x = 1. */
        const geomlib_vec3 p[6] = {
            {0,0,0},{1,0,0},{1,1,0},{0,1,0},{1,0,-1},{1,1,-1}
        };
        for (int i = 0; i < 6; ++i)
            geomlib_trimesh_add_vertex(&bent, p[i], nullptr, nullptr);
        geomlib_trimesh_add_triangle(&bent, 0, 1, 2);
        geomlib_trimesh_add_triangle(&bent, 0, 2, 3);
        geomlib_trimesh_add_triangle(&bent, 1, 4, 5);
        geomlib_trimesh_add_triangle(&bent, 1, 5, 2);

        const geomlib_atlas_segmentation segs[] = {
            GEOMLIB_ATLAS_SEG_PLANARITY, GEOMLIB_ATLAS_SEG_DISTORTION
        };
        const geomlib_atlas_packing packs[] = {
            GEOMLIB_ATLAS_PACK_SHELF, GEOMLIB_ATLAS_PACK_MAXRECTS,
            GEOMLIB_ATLAS_PACK_GUILLOTINE
        };
        for (int s = 0; s < 2; ++s) {
            for (int k = 0; k < 3; ++k) {
                geomlib_atlas_params params = geomlib_atlas_params_default();
                params.segmentation = segs[s];
                params.packing = packs[k];
                params.parameterization = GEOMLIB_PARAM_HARMONIC;
                params.optimize_rotation = true;
                geomlib_atlas a = {};
                a.charts = (geomlib_atlas_chart*)calloc(256, sizeof(geomlib_atlas_chart));
                CHECK(geomlib_generate_atlas(&bent, &params, &a) == GEOMLIB_MESH_SUCCESS,
                      "atlas pipeline succeeds");
                CHECK(a.chart_count >= 1 && a.chart_count <= 256, "chart count in range");
                CHECK(a.coverage > 0.0f && a.coverage <= 100.0f, "coverage is a percentage");
                geomlib_atlas_free(&a);
            }
        }
        geomlib_trimesh_free(&bent);
    }

    std::printf("\n=== ALL TESTS PASSED ===\n");
    return 0;
}
