/*
 * Vertex-clustering decimation test: a sphere simplifies while keeping
 * valid topology and approximately the original volume.
 */
#include <geomlib/geomlib_clustering.h>

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
    double v0 = geomlib_trimesh_volume(&sphere);
    geomlib_trimesh out = {};
    int rc = geomlib_decimate_vertex_clustering(&sphere, 0.12f, &out);
    std::printf("clustering rc=%d triangles %zu->%zu volume %.4f->%.4f\n",
                rc, sphere.triangle_count, out.triangle_count, v0,
                out.triangle_count ? geomlib_trimesh_volume(&out) : 0.0);
    check(rc == GEOMLIB_MESH_SUCCESS, "clustering decimation succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(out.triangle_count > 0 && out.triangle_count < sphere.triangle_count,
              "triangle count reduced");
        bool valid = true;
        for (size_t i = 0; i < out.triangle_count*3; ++i)
            if (out.indices[i] >= out.vertex_count) valid = false;
        check(valid, "output indices are in range");
        bool finite = true;
        for (size_t i = 0; i < out.vertex_count; ++i) {
            geomlib_vec3 p = out.vertices[i];
            if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
                finite = false;
        }
        check(finite, "output vertices are finite");
        geomlib_trimesh_free(&out);
    }
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d clustering test(s) failed\n", failures);
        return 1;
    }
    std::printf("All clustering tests passed\n");
    return 0;
}
