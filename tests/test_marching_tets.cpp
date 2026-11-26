/*
 * Marching tetrahedra test: extract a unit sphere isosurface and check
 * volume and watertightness.
 */
#include <geomlib/geomlib_marching_tets.h>
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
    geomlib_sdf sdf = geomlib_mesh_to_sdf_resolution(&sphere, 32);
    check(sdf.values != NULL, "SDF generated");
    if (sdf.values) {
        geomlib_trimesh out = {};
        int rc = geomlib_sdf_marching_tets(&sdf, &out);
        double vol = rc == GEOMLIB_MESH_SUCCESS ? geomlib_trimesh_volume(&out) : 0.0;
        geomlib_manifold_analysis a = {};
        if (rc == GEOMLIB_MESH_SUCCESS) a = geomlib_analyze_manifold(&out);
        std::printf("marching tets rc=%d triangles=%zu volume=%.4f closed=%d\n",
                    rc, out.triangle_count, vol, a.is_closed);
        check(rc == GEOMLIB_MESH_SUCCESS && out.triangle_count > 100,
              "marching tetrahedra produces a mesh");
        check(std::fabs(fabs(vol) - 4.18879) < 0.5, "sphere volume is approximately correct");
        check(a.is_closed, "extracted surface is watertight");
        geomlib_manifold_analysis_free(&a);
        geomlib_trimesh_free(&out);
    }
    geomlib_sdf_free(&sdf);
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d marching tetrahedra test(s) failed\n", failures);
        return 1;
    }
    std::printf("All marching tetrahedra tests passed\n");
    return 0;
}
