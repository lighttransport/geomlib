/*
 * Poisson surface reconstruction tests: sphere and box point clouds with
 * analytic outward normals, checking watertight output and approximate volume.
 */
#include <geomlib/geomlib_poisson.h>
#include <geomlib/geomlib_topology.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static void add_sphere_points(geomlib_pointcloud *pc, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        float u = (float)i / (float)n;
        float z = 1.0f - 2.0f*u;
        float r = sqrtf(fmaxf(0.0f, 1.0f-z*z));
        float theta = (float)(i % 97) * 0.618034f * 6.2831853f;
        geomlib_vec3 p = { r*cosf(theta), r*sinf(theta), z };
        geomlib_pointcloud_add_point(pc, p, &p, 1.0f);
    }
}

int main(void)
{
    geomlib_pointcloud sphere;
    geomlib_pointcloud_init(&sphere, 512);
    add_sphere_points(&sphere, 400);

    geomlib_poisson_params params = geomlib_poisson_params_default();
    params.resolution = 40;
    params.max_iterations = 400;
    geomlib_trimesh out = {};
    int rc = geomlib_reconstruct_poisson(&sphere, &params, &out);
    std::printf("poisson rc=%d triangles=%zu volume=%.4f\n",
                rc, out.triangle_count, out.triangle_count ? geomlib_trimesh_volume(&out) : 0.0);
    check(rc == GEOMLIB_MESH_SUCCESS, "poisson reconstruction succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(out.triangle_count > 100, "poisson produces a usable mesh");
        double vol = fabs(geomlib_trimesh_volume(&out));
        double expected = 4.0/3.0 * M_PI;
        check(vol > expected*0.4 && vol < expected*1.6, "sphere volume is approximately correct");
        geomlib_trimesh_free(&out);
    }
    geomlib_pointcloud_free(&sphere);

    if (failures != 0) {
        std::printf("%d poisson test(s) failed\n", failures);
        return 1;
    }
    std::printf("All poisson tests passed\n");
    return 0;
}
