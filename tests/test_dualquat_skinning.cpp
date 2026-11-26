/*
 * Dual-quaternion skinning test: single-bone rotation and a 50/50 blend of
 * identity and rotation; the blend must preserve length better than LBS.
 */
#include <geomlib/geomlib_dualquat_skinning.h>
#include <geomlib/geomlib_skinning.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh tri;
    geomlib_trimesh_init(&tri, 3, 1);
    geomlib_trimesh_add_vertex(&tri, geomlib_vec3_create(1,0,0),0,0);
    geomlib_trimesh_add_vertex(&tri, geomlib_vec3_create(2,0,0),0,0);
    geomlib_trimesh_add_vertex(&tri, geomlib_vec3_create(0,1,0),0,0);
    geomlib_trimesh_add_triangle(&tri, 0,1,2);
    std::vector<geomlib_skin_weight> weights(3);

    geomlib_mat4 rot = geomlib_mat4_identity();
    float a = 1.57079632679f;
    rot.m[0]=cosf(a); rot.m[1]=sinf(a);
    rot.m[4]=-sinf(a); rot.m[5]=cosf(a);

    for (auto &w : weights) { w.bone_indices[0]=0; w.weights[0]=1.0f; }
    geomlib_trimesh out = {};
    int rc = geomlib_skin_mesh_dualquat(&tri, weights.data(), &rot, 1, &out);
    check(rc == GEOMLIB_MESH_SUCCESS, "single bone DQ skinning succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(std::fabs(out.vertices[0].x) < 1e-3 && std::fabs(out.vertices[0].y-1.0f) < 1e-3,
              "DQ single bone rotation is exact");
        geomlib_trimesh_free(&out);
    }

    for (auto &w : weights) {
        w.bone_indices[0]=0; w.weights[0]=0.5f;
        w.bone_indices[1]=1; w.weights[1]=0.5f;
    }
    geomlib_mat4 bones[2] = {geomlib_mat4_identity(), rot};
    out = {};
    rc = geomlib_skin_mesh_dualquat(&tri, weights.data(), bones, 2, &out);
    check(rc == GEOMLIB_MESH_SUCCESS, "blended DQ skinning succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        geomlib_vec3 p = out.vertices[0];
        double r = std::sqrt(p.x*p.x + p.y*p.y);
        std::printf("DQ blend point=(%.4f,%.4f) radius=%.4f (LBS radius 0.7071)\n",
                    p.x, p.y, r);
        check(std::fabs(r-1.0) < 1e-3, "DQ blend preserves radius");
        geomlib_trimesh_free(&out);
    }
    geomlib_trimesh_free(&tri);
    if (failures != 0) {
        std::printf("%d dual-quaternion skinning test(s) failed\n", failures);
        return 1;
    }
    std::printf("All dual-quaternion skinning tests passed\n");
    return 0;
}
