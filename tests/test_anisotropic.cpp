/*
 * Anisotropic remeshing tests: a sphere should gain vertices (refinement),
 * stay finite, and remain near the original surface.
 */
#include <geomlib/geomlib_anisotropic.h>

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
    geomlib_create_sphere(&sphere, 1.0f, 16, 8);
    size_t v0 = sphere.vertex_count, t0 = sphere.triangle_count;

    geomlib_anisotropic_params params = geomlib_anisotropic_params_default(0.18f);
    params.iterations = 3;
    int rc = geomlib_remesh_anisotropic(&sphere, &params);
    std::printf("anisotropic rc=%d vertices %zu->%zu triangles %zu->%zu\n",
                rc, v0, sphere.vertex_count, t0, sphere.triangle_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "anisotropic remeshing succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        bool finite = true, on_surface = true;
        for (size_t i = 0; i < sphere.vertex_count; ++i) {
            geomlib_vec3 p = sphere.vertices[i];
            if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
                finite = false;
            float r = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
            if (std::fabs(r - 1.0f) > 0.1f) on_surface = false;
        }
        check(finite, "all vertices are finite");
        check(on_surface, "vertices remain near the unit sphere");
        check(sphere.vertex_count >= v0, "refinement does not lose vertices");
    }
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d anisotropic test(s) failed\n", failures);
        return 1;
    }
    std::printf("All anisotropic tests passed\n");
    return 0;
}
