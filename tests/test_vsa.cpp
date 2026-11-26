/*
 * Variational shape approximation tests: sphere clustering and simplified
 * output mesh.
 */
#include <geomlib/geomlib_vsa.h>
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
    geomlib_create_sphere(&sphere, 1.0f, 24, 12);
    geomlib_vsa_result result = geomlib_vsa_result_create();
    geomlib_trimesh out = {};
    int rc = geomlib_vsa_simplify(&sphere, 12, 5, &result, &out);
    std::printf("vsa rc=%d regions=%zu triangles %zu->%zu\n",
                rc, result.region_count, sphere.triangle_count, out.triangle_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "VSA succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(result.region_count == 12 && result.triangle_count == sphere.triangle_count,
              "clustering dimensions");
        bool valid = true;
        for (size_t i = 0; i < result.triangle_count; ++i)
            if (result.region_of_triangle[i] >= result.region_count) valid = false;
        check(valid, "region ids are in range");
        check(out.triangle_count > 0 && out.triangle_count < sphere.triangle_count,
              "output mesh is simplified");
        geomlib_trimesh_free(&out);
    }
    geomlib_vsa_result_free(&result);
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d VSA test(s) failed\n", failures);
        return 1;
    }
    std::printf("All VSA tests passed\n");
    return 0;
}
