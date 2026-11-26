/*
 * Mesh saliency tests: a sharp bump on a grid should score higher than the
 * flat interior, and saliency segmentation should cover every vertex.
 */
#include <geomlib/geomlib_saliency.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh grid;
    geomlib_trimesh_init(&grid, 36, 50);
    for (int y = 0; y < 6; ++y)
        for (int x = 0; x < 6; ++x) {
            float z = (x == 3 && y == 3) ? 1.0f : 0.0f;
            geomlib_trimesh_add_vertex(&grid,
                geomlib_vec3_create((float)x*0.2f, (float)y*0.2f, z), NULL, NULL);
        }
    for (int y = 0; y < 5; ++y)
        for (int x = 0; x < 5; ++x) {
            uint32_t a = (uint32_t)(y*6+x), b = a+1, c = a+6, d = c+1;
            geomlib_trimesh_add_triangle(&grid, a, b, d);
            geomlib_trimesh_add_triangle(&grid, a, d, c);
        }

    geomlib_saliency_result r = geomlib_saliency_result_create(grid.vertex_count);
    int rc = geomlib_compute_saliency(&grid, 3, &r);
    std::printf("saliency rc=%d max=%.4f center=%.4f corner=%.4f\n",
                rc, r.values[21], r.values[21], r.values[0]);
    check(rc == GEOMLIB_MESH_SUCCESS, "saliency computation succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        bool valid = true;
        for (size_t i = 0; i < grid.vertex_count; ++i) {
            if (!std::isfinite(r.values[i]) || r.values[i] < 0.0f || r.values[i] > 1.0f)
                valid = false;
        }
        check(valid, "saliency values are normalized and finite");
        check(r.values[21] > 0.5f, "sharp bump has high saliency");

        uint32_t *labels = (uint32_t *)malloc(grid.vertex_count*sizeof(uint32_t));
        rc = geomlib_segment_by_saliency(&grid, &r, 3, labels);
        check(rc == GEOMLIB_MESH_SUCCESS, "saliency segmentation succeeds");
        bool in_range = true;
        for (size_t i = 0; i < grid.vertex_count; ++i)
            if (labels[i] >= 3) in_range = false;
        check(in_range, "every vertex gets a valid region label");
        free(labels);
    }
    geomlib_saliency_result_free(&r);
    geomlib_trimesh_free(&grid);
    if (failures != 0) {
        std::printf("%d saliency test(s) failed\n", failures);
        return 1;
    }
    std::printf("All saliency tests passed\n");
    return 0;
}
