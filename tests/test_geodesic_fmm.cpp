/*
 * Fast-marching geodesic test: distances on a unit sphere from one vertex.
 */
#include <geomlib/geomlib_geodesic_fmm.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh sphere;
    geomlib_create_sphere(&sphere, 1.0f, 32, 16);
    geomlib_geodesic_fmm_result r = geomlib_geodesic_fmm_result_create(0);
    int rc = geomlib_geodesic_fast_marching(&sphere, 0, 6, &r);
    float mx = 0.0f;
    for (size_t i = 0; i < r.count; ++i) mx = std::max(mx, r.distances[i]);
    std::printf("geodesic fmm rc=%d count=%zu source=%.4f max=%.4f\n",
                rc, r.count, r.distances[0], mx);
    check(rc == GEOMLIB_MESH_SUCCESS, "fast marching succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(r.count == sphere.vertex_count, "one distance per vertex");
        check(std::fabs(r.distances[0]) < 1e-6f, "source distance is zero");
        bool finite = true;
        for (size_t i = 0; i < r.count; ++i)
            if (!std::isfinite(r.distances[i])) finite = false;
        check(finite, "all distances are finite");
        check(mx > 2.5f && mx < 4.0f, "sphere diameter is near pi");
    }
    geomlib_geodesic_fmm_result_free(&r);
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d geodesic FMM test(s) failed\n", failures);
        return 1;
    }
    std::printf("All geodesic FMM tests passed\n");
    return 0;
}
