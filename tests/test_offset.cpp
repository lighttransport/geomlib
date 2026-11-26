/*
 * Mesh offset tests: outward/inward offsets of a sphere change its radius and
 * enclosed volume in the expected direction.
 */
#include <geomlib/geomlib_offset.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static float radius_of(const geomlib_trimesh &m)
{
    geomlib_aabb b = geomlib_trimesh_compute_aabb(&m);
    return 0.5f*(geomlib_aabb_size(b).x + geomlib_aabb_size(b).y +
                 geomlib_aabb_size(b).z)/3.0f;
}

int main(void)
{
    geomlib_trimesh sphere = {};
    int rc = geomlib_create_sphere(&sphere, 0.5f, 24, 16);
    check(rc == GEOMLIB_MESH_SUCCESS, "sphere created");
    if (rc != GEOMLIB_MESH_SUCCESS) return 1;

    geomlib_offset_params p = geomlib_offset_params_default();
    p.resolution = 40;
    p.padding = 6.0f;

    geomlib_trimesh out = {};
    rc = geomlib_offset_mesh(&sphere, 0.1f, &p, &out);
    std::printf("offset +0.1 rc=%d tris=%zu radius=%.4f volume=%.4f\n",
                rc, out.triangle_count, rc == GEOMLIB_MESH_SUCCESS ? radius_of(out) : 0.0f,
                rc == GEOMLIB_MESH_SUCCESS ? geomlib_trimesh_volume(&out) : 0.0f);
    check(rc == GEOMLIB_MESH_SUCCESS && out.triangle_count > 0, "outward offset builds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        float r = radius_of(out);
        float v = geomlib_trimesh_volume(&out);
        check(r > 0.55f && r < 0.65f, "outward radius grows to about 0.6");
        check(v > 0.7f && v < 1.1f, "outward volume near analytic 0.905");
        check(v > geomlib_trimesh_volume(&sphere), "outward volume increases");
        geomlib_trimesh_free(&out);
    }

    rc = geomlib_offset_mesh(&sphere, -0.1f, &p, &out);
    std::printf("offset -0.1 rc=%d tris=%zu radius=%.4f volume=%.4f\n",
                rc, out.triangle_count, rc == GEOMLIB_MESH_SUCCESS ? radius_of(out) : 0.0f,
                rc == GEOMLIB_MESH_SUCCESS ? geomlib_trimesh_volume(&out) : 0.0f);
    check(rc == GEOMLIB_MESH_SUCCESS && out.triangle_count > 0, "inward offset builds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        float r = radius_of(out);
        float v = geomlib_trimesh_volume(&out);
        check(r > 0.35f && r < 0.45f, "inward radius shrinks to about 0.4");
        check(v > 0.15f && v < 0.40f, "inward volume near analytic 0.268");
        geomlib_trimesh_free(&out);
    }

    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d offset test(s) failed\n", failures);
        return 1;
    }
    std::printf("All offset tests passed\n");
    return 0;
}
