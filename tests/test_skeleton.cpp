/*
 * Curve skeleton tests: a sphere collapses to a small graph with valid
 * node/edge indices and finite positions; the input mesh is unchanged.
 */
#include <geomlib/geomlib_skeleton.h>

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
    geomlib_create_sphere(&sphere, 1.0f, 20, 10);
    size_t v0 = sphere.vertex_count, t0 = sphere.triangle_count;

    geomlib_skeleton_params params = geomlib_skeleton_params_default();
    params.collapse_ratio = 0.05f;
    params.max_iterations = 100;
    geomlib_curve_skeleton skeleton = geomlib_curve_skeleton_create();
    int rc = geomlib_extract_curve_skeleton(&sphere, &params, &skeleton);
    std::printf("skeleton rc=%d nodes=%zu edges=%zu\n",
                rc, skeleton.point_count, skeleton.edge_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "skeleton extraction succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(skeleton.point_count > 0 && skeleton.point_count < v0,
              "skeleton is smaller than the input");
        check(skeleton.edge_count > 0, "skeleton has edges");
        bool valid = true;
        for (size_t i = 0; i < skeleton.point_count; ++i) {
            geomlib_vec3 p = skeleton.points[i];
            if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
                valid = false;
            if (skeleton.original_vertices[i] >= v0) valid = false;
        }
        for (size_t i = 0; i < skeleton.edge_count*2; ++i)
            if (skeleton.edges[i] >= skeleton.point_count) valid = false;
        check(valid, "skeleton nodes and edges are valid");
    }
    check(sphere.vertex_count == v0 && sphere.triangle_count == t0,
          "input mesh is unchanged");
    geomlib_curve_skeleton_free(&skeleton);
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d skeleton test(s) failed\n", failures);
        return 1;
    }
    std::printf("All skeleton tests passed\n");
    return 0;
}
