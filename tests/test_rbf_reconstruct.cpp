/*
 * RBF reconstruction test: sphere point cloud yields a closed surface.
 */
#include <geomlib/geomlib_rbf_reconstruct.h>
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
    geomlib_pointcloud cloud;
    geomlib_pointcloud_init(&cloud, 256);
    for (int i = 0; i < 200; ++i) {
        float u = i*0.618033f*6.283185f;
        float z = 1.0f - 2.0f*(float)(i+0.5)/200.0f;
        float r = sqrtf(fmaxf(0.0f, 1.0f-z*z));
        geomlib_vec3 p = geomlib_vec3_create(r*cosf(u), r*sinf(u), z);
        geomlib_pointcloud_add_point(&cloud, p, &p, 1.0f);
    }
    geomlib_rbf_params params = geomlib_rbf_params_default();
    params.resolution = 32;
    params.max_points = 200;
    geomlib_trimesh out = {};
    int rc = geomlib_reconstruct_rbf(&cloud, &params, &out);
    double vol = rc == GEOMLIB_MESH_SUCCESS ? geomlib_trimesh_volume(&out) : 0.0;
    std::printf("rbf rc=%d triangles=%zu volume=%.4f (expected ~4.19)\n",
                rc, out.triangle_count, vol);
    check(rc == GEOMLIB_MESH_SUCCESS, "RBF reconstruction succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(out.triangle_count > 100, "surface has triangles");
        check(std::fabs(fabs(vol)-4.18879) < 1.5, "sphere volume is approximately correct");
        geomlib_trimesh_free(&out);
    }
    geomlib_pointcloud_free(&cloud);
    if (failures != 0) {
        std::printf("%d RBF reconstruction test(s) failed\n", failures);
        return 1;
    }
    std::printf("All RBF reconstruction tests passed\n");
    return 0;
}
