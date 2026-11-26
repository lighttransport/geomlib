/*
 * MLS projection tests: planar and spherical point clouds.
 */
#include <geomlib/geomlib_mls.h>
#include <geomlib/geomlib_kdtree.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_pointcloud plane;
    geomlib_pointcloud_init(&plane, 400);
    for (int i = 0; i < 400; ++i) {
        float x = ((i*37)%101)/100.0f*2.0f-1.0f;
        float y = ((i*73)%97)/100.0f*2.0f-1.0f;
        geomlib_pointcloud_add_point(&plane, geomlib_vec3_create(x,y,0),
                                     nullptr, 1.0f);
    }
    geomlib_mls_params p = geomlib_mls_params_default();
    p.neighbors = 24;
    geomlib_vec3 out, nrm;
    int rc = geomlib_mls_project(&plane, geomlib_vec3_create(0.1f,0.2f,0.3f),
                                 &p, &out, &nrm);
    std::printf("mls plane rc=%d z=%.4f normal=(%.2f,%.2f,%.2f)\n",
                rc, out.z, nrm.x, nrm.y, nrm.z);
    check(rc == GEOMLIB_MESH_SUCCESS, "plane projection succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(std::fabs(out.z) < 0.02f, "plane projection lands on z=0");
        check(std::fabs(std::fabs(nrm.z)-1.0f) < 0.02f, "plane normal is vertical");
    }
    geomlib_pointcloud_free(&plane);

    geomlib_pointcloud sphere;
    geomlib_pointcloud_init(&sphere, 500);
    for (int i = 0; i < 500; ++i) {
        float u = i*0.618033f*6.283185f;
        float z = 1.0f-2.0f*(float)(i+0.5)/500.0f;
        float r = sqrtf(fmaxf(0.0f,1.0f-z*z));
        geomlib_pointcloud_add_point(&sphere,
            geomlib_vec3_create(r*cosf(u),r*sinf(u),z), nullptr, 1.0f);
    }
    geomlib_vec3 q = geomlib_vec3_create(1.2f,0.0f,0.0f);
    rc = geomlib_mls_project(&sphere, q, &p, &out, &nrm);
    float rad = geomlib_vec3_length(out);
    std::printf("mls sphere rc=%d radius=%.4f (query %.4f)\n",
                rc, rad, geomlib_vec3_length(q));
    check(rc == GEOMLIB_MESH_SUCCESS, "sphere projection succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        geomlib_kdtree *tree = geomlib_kdtree_build(sphere.points, sphere.point_count);
        uint32_t id=0;
        float db=0.0f, da=0.0f;
        geomlib_kdtree_knn(tree, q, 1, &id, &db);
        geomlib_kdtree_knn(tree, out, 1, &id, &da);
        geomlib_kdtree_free(tree);
        std::printf("mls sphere distance before=%.4f after=%.4f\n", db, da);
        check(da < db, "projection moves toward the point cloud");
        check(rad > 0.9f && rad < 1.1f, "projected radius is near 1");
        check(std::fabs(geomlib_vec3_length(nrm)-1.0f) < 1e-3f,
              "sphere normal is unit length");
    }
    geomlib_pointcloud_free(&sphere);

    if (failures != 0) {
        std::printf("%d MLS test(s) failed\n", failures);
        return 1;
    }
    std::printf("All MLS tests passed\n");
    return 0;
}
