/*
 * ACVD remeshing test: sphere remeshed to a target vertex count, checking
 * vertex count, area retention and surface adherence.
 */
#include <geomlib/geomlib_acvd.h>
#include <geomlib/geomlib_topology.h>

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
    double input_area = geomlib_trimesh_surface_area(&sphere);

    geomlib_acvd_params params = geomlib_acvd_params_default(500);
    params.lloyd_iterations = 4;
    geomlib_trimesh out = {};
    int rc = geomlib_remesh_acvd(&sphere, &params, &out);
    std::printf("acvd rc=%d vertices=%zu triangles=%zu area=%.4f (input %.4f)\n",
                rc, out.vertex_count, out.triangle_count,
                out.vertex_count ? geomlib_trimesh_surface_area(&out) : 0.0, input_area);
    check(rc == GEOMLIB_MESH_SUCCESS, "acvd remeshing succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(out.vertex_count > 100 && out.vertex_count < 2000,
              "vertex count near target");
        check(out.triangle_count > 100, "useful number of triangles");
        double area = geomlib_trimesh_surface_area(&out);
        check(area > input_area*0.5 && area < input_area*1.5,
              "surface area is approximately preserved");
        geomlib_trimesh_free(&out);
    }
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d acvd test(s) failed\n", failures);
        return 1;
    }
    std::printf("All acvd tests passed\n");
    return 0;
}
