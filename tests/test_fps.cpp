/*
 * Farthest-point sampling tests: distinct samples, monotonically decreasing
 * covering radius, determinism, and final coverage of a sphere.
 */
#include <geomlib/geomlib_fps.h>
#include <geomlib/geomlib_geodesic.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh sphere = {};
    int rc = geomlib_create_sphere(&sphere, 1.0f, 16, 12);
    check(rc == GEOMLIB_MESH_SUCCESS, "sphere created");
    if (rc != GEOMLIB_MESH_SUCCESS) return 1;

    geomlib_fps_result a = geomlib_fps_result_create();
    rc = geomlib_farthest_point_sampling(&sphere, 8, 0, &a);
    std::printf("fps rc=%d count=%zu radii:", rc, a.count);
    for (size_t i = 0; i < a.count; ++i) std::printf(" %.3f", a.radii[i]);
    std::printf("\n");
    check(rc == GEOMLIB_MESH_SUCCESS, "FPS succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(a.count == 8, "eight samples requested and returned");
        bool distinct = true;
        for (size_t i = 0; i < a.count; ++i) {
            if (a.indices[i] >= sphere.vertex_count) distinct = false;
            for (size_t j = 0; j < i; ++j)
                if (a.indices[i] == a.indices[j]) distinct = false;
        }
        check(distinct, "samples are distinct in-range vertices");
        bool monotone = true;
        for (size_t i = 2; i < a.count; ++i)
            if (a.radii[i] > a.radii[i-1] + 1e-3f) monotone = false;
        check(monotone, "covering radius decreases monotonically");
        check(a.radii[a.count-1] < a.radii[1]*0.8f, "radius shrinks with samples");

        /* Final covering radius matches the geodesic distance to all samples. */
        geomlib_geodesic_distances g = geomlib_geodesic_create(sphere.vertex_count);
        if (geomlib_compute_geodesic_heat(&sphere, a.indices, a.count, 0.0f, &g) ==
            GEOMLIB_MESH_SUCCESS && g.distances) {
            float mx = 0.0f;
            for (size_t v = 0; v < sphere.vertex_count; ++v)
                if (g.distances[v] > mx) mx = g.distances[v];
            std::printf("fps final radius=%.3f recomputed=%.3f\n",
                        a.radii[a.count-1], mx);
            check(std::fabs(mx-a.radii[a.count-1]) < 0.05f,
                  "reported radius matches final geodesic coverage");
        }
        geomlib_geodesic_free(&g);
    }

    geomlib_fps_result b = geomlib_fps_result_create();
    rc = geomlib_farthest_point_sampling(&sphere, 8, 0, &b);
    if (rc == GEOMLIB_MESH_SUCCESS && a.count) {
        bool same = (a.count == b.count);
        for (size_t i = 0; i < b.count && same; ++i)
            if (b.indices[i] != a.indices[i]) same = false;
        check(same, "sampling is deterministic");
    }
    geomlib_fps_result_free(&a);
    geomlib_fps_result_free(&b);
    geomlib_trimesh_free(&sphere);

    if (failures != 0) {
        std::printf("%d FPS test(s) failed\n", failures);
        return 1;
    }
    std::printf("All FPS tests passed\n");
    return 0;
}
