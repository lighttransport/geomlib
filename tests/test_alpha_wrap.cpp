/*
 * Alpha wrapping tests: an open surface becomes a watertight enclosing mesh.
 */
#include <geomlib/geomlib_alpha_wrap.h>
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
    geomlib_trimesh plane;
    geomlib_create_plane(&plane, 1.0f, 1.0f, 4);
    geomlib_alpha_wrap_params params = geomlib_alpha_wrap_params_default();
    params.resolution = 32;
    params.alpha = 0.08f;
    params.smoothing_steps = 1;
    geomlib_trimesh out = {};
    int rc = geomlib_alpha_wrap(&plane, &params, &out);
    std::printf("alpha wrap rc=%d triangles=%zu volume=%.4f\n",
                rc, out.triangle_count,
                out.triangle_count ? geomlib_trimesh_volume(&out) : 0.0);
    check(rc == GEOMLIB_MESH_SUCCESS, "alpha wrap succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        geomlib_manifold_analysis a = geomlib_analyze_manifold(&out);
        check(out.triangle_count > 0, "wrap has triangles");
        if (!(a.is_closed && a.is_manifold))
            std::printf("  closed=%d manifold=%d edge=%d vertex=%d bad=%zu\n",
                        a.is_closed?1:0, a.is_manifold?1:0,
                        a.is_edge_manifold?1:0, a.is_vertex_manifold?1:0,
                        a.nonmanifold_count);
        check(a.is_closed && a.is_manifold, "wrap is watertight and manifold");
        check(std::fabs(geomlib_trimesh_volume(&out)) > 0.0f, "wrap encloses a volume");
        geomlib_manifold_analysis_free(&a);
        geomlib_trimesh_free(&out);
    }
    geomlib_trimesh_free(&plane);
    if (failures != 0) {
        std::printf("%d alpha wrap test(s) failed\n", failures);
        return 1;
    }
    std::printf("All alpha wrap tests passed\n");
    return 0;
}
