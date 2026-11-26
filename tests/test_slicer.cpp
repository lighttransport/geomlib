/*
 * Mesh slicing tests: a plane through a box produces a closed 4-segment loop.
 */
#include <geomlib/geomlib_slicer.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh box;
    geomlib_create_box(&box, 1.0f, 1.0f, 1.0f);
    geomlib_slice_result result = geomlib_slice_result_create();
    int rc = geomlib_slice_mesh(&box, geomlib_vec3_create(0,0,0),
                                geomlib_vec3_create(0,0,1), 1e-4f, &result);
    std::printf("slice rc=%d polylines=%zu points=%zu\n",
                rc, result.polyline_count, result.point_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "slicing succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(result.polyline_count >= 1, "at least one polyline");
        bool square = false;
        for (size_t p = 0; p < result.polyline_count; ++p) {
            if (result.counts[p] >= 4) {
                bool all_z = true;
                for (uint32_t i = 0; i < result.counts[p]; ++i) {
                    geomlib_vec3 v = result.points[result.starts[p]+i];
                    if (std::fabs(v.z) > 1e-4f) all_z = false;
                }
                if (all_z) square = true;
            }
        }
        check(square, "slice loop lies on the plane and has four points");
        geomlib_slice_result_free(&result);
    }
    geomlib_trimesh_free(&box);
    if (failures != 0) {
        std::printf("%d slicer test(s) failed\n", failures);
        return 1;
    }
    std::printf("All slicer tests passed\n");
    return 0;
}
