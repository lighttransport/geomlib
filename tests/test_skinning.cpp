/*
 * Linear blend skinning test: a 90-degree bone rotation moves weighted
 * vertices exactly, and blended weights produce the midpoint.
 */
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
    geomlib_trimesh_add_vertex(&tri, geomlib_vec3_create(1,0,0), 0, 0);
    geomlib_trimesh_add_vertex(&tri, geomlib_vec3_create(2,0,0), 0, 0);
    geomlib_trimesh_add_vertex(&tri, geomlib_vec3_create(0,1,0), 0, 0);
    geomlib_trimesh_add_triangle(&tri, 0, 1, 2);

    std::vector<geomlib_skin_weight> weights(3);
    for (auto &w : weights) { w.bone_indices[0]=0; w.weights[0]=1.0f; }

    geomlib_mat4 rot = geomlib_mat4_identity();
    float a = 1.57079632679f;
    rot.m[0] = cosf(a); rot.m[1] = sinf(a);
    rot.m[4] = -sinf(a); rot.m[5] = cosf(a);

    geomlib_trimesh out = {};
    int rc = geomlib_skin_mesh(&tri, weights.data(), &rot, 1, &out);
    std::printf("skinning rc=%d v0=(%.4f,%.4f)\n", rc,
                out.vertices[0].x, out.vertices[0].y);
    check(rc == GEOMLIB_MESH_SUCCESS, "skinning succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(std::fabs(out.vertices[0].x) < 1e-4 && std::fabs(out.vertices[0].y-1.0f) < 1e-4,
              "single-bone rotation is exact");
        check(std::fabs(out.vertices[2].x+1.0f) < 1e-4, "second vertex rotates");
        geomlib_trimesh_free(&out);
    }

    /* Blend two transforms: identity and 90-degree rotation, weights 0.5. */
    for (auto &w : weights) {
        w.bone_indices[0]=0; w.weights[0]=0.5f;
        w.bone_indices[1]=1; w.weights[1]=0.5f;
    }
    geomlib_mat4 bones[2] = {geomlib_mat4_identity(), rot};
    out = {};
    rc = geomlib_skin_mesh(&tri, weights.data(), bones, 2, &out);
    check(rc == GEOMLIB_MESH_SUCCESS, "blended skinning succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        geomlib_vec3 p = out.vertices[0];
        check(std::fabs(p.x-0.5f) < 1e-4 && std::fabs(p.y-0.5f) < 1e-4,
              "blend produces the midpoint");
        geomlib_trimesh_free(&out);
    }
    geomlib_trimesh_free(&tri);
    if (failures != 0) {
        std::printf("%d skinning test(s) failed\n", failures);
        return 1;
    }
    std::printf("All skinning tests passed\n");
    return 0;
}
