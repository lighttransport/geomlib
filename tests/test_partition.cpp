/*
 * Mesh partition tests: a sphere splits into four balanced parts covering
 * every vertex.
 */
#include <geomlib/geomlib_partition.h>

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
    geomlib_partition_result r = geomlib_partition_result_create();
    int rc = geomlib_partition_mesh(&sphere, 4, &r);
    std::printf("partition rc=%d parts=%zu sizes:", rc, r.part_count);
    for (size_t i = 0; i < r.part_count; ++i) std::printf(" %zu", r.part_sizes[i]);
    std::printf("\n");
    check(rc == GEOMLIB_MESH_SUCCESS, "partition succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        size_t total = 0;
        bool valid = true;
        for (size_t i = 0; i < r.vertex_count; ++i)
            if (r.part_of_vertex[i] >= r.part_count) valid = false;
        check(valid, "all labels are in range");
        for (size_t i = 0; i < r.part_count; ++i) total += r.part_sizes[i];
        check(total == r.vertex_count, "part sizes sum to vertex count");
        bool balanced = true;
        for (size_t i = 0; i < r.part_count; ++i)
            if (r.part_sizes[i] > (r.vertex_count+3)/4 + 4) balanced = false;
        check(balanced, "parts are approximately balanced");
        geomlib_partition_result_free(&r);
    }
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d partition test(s) failed\n", failures);
        return 1;
    }
    std::printf("All partition tests passed\n");
    return 0;
}
