/*
 * Convex decomposition tests: a box decomposes into convex pieces with
 * positive volume.
 */
#include <geomlib/geomlib_convex_decomposition.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh box;
    geomlib_create_box(&box, 1.0f, 1.0f, 1.0f);
    geomlib_convex_decomposition d = geomlib_convex_decomposition_create();
    int rc = geomlib_approx_convex_decomposition(&box, 4, &d);
    std::printf("convex decomposition rc=%d pieces=%zu\n", rc, d.piece_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "decomposition succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(d.piece_count == 4, "requested piece count");
        int positive = 0;
        for (size_t i = 0; i < d.piece_count; ++i) {
            double v = geomlib_trimesh_volume(&d.pieces[i]);
            if (v > 0.0 && d.pieces[i].triangle_count > 0) positive++;
        }
        check(positive >= 2, "at least two convex pieces have positive volume");
        geomlib_convex_decomposition_free(&d);
    }
    geomlib_trimesh_free(&box);
    if (failures != 0) {
        std::printf("%d convex decomposition test(s) failed\n", failures);
        return 1;
    }
    std::printf("All convex decomposition tests passed\n");
    return 0;
}
