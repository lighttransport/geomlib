/*
 * Bounded Hausdorff tests: bound must contain the sampled value, vanish for
 * identical meshes, and grow when a mesh is translated.
 */
#include <geomlib/geomlib_distance.h>
#include <geomlib/geomlib_mesh.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static geomlib_trimesh make_box(float sx)
{
    geomlib_trimesh m;
    geomlib_trimesh_init(&m, 8, 12);
    float v[8][3] = {
        {0,0,0},{sx,0,0},{sx,sx,0},{0,sx,0},
        {0,0,sx},{sx,0,sx},{sx,sx,sx},{0,sx,sx}
    };
    int f[12][3] = {
        {0,1,2},{0,2,3},{4,6,5},{4,7,6},
        {0,4,5},{0,5,1},{3,2,6},{3,6,7},
        {0,3,7},{0,7,4},{1,5,6},{1,6,2}
    };
    for (int i = 0; i < 8; ++i)
        geomlib_trimesh_add_vertex(&m, geomlib_vec3_create(v[i][0],v[i][1],v[i][2]), 0, 0);
    for (int i = 0; i < 12; ++i)
        geomlib_trimesh_add_triangle(&m, f[i][0], f[i][1], f[i][2]);
    return m;
}

int main(void)
{
    geomlib_trimesh a = make_box(1.0f);
    geomlib_trimesh b = make_box(1.0f);
    geomlib_hausdorff_result r;

    check(geomlib_hausdorff_bounded(&a, &b, &r) == GEOMLIB_MESH_SUCCESS, "bounded rc");
    check(r.samples > 0, "samples counted");
    check(r.lower_bound <= r.upper_bound + 1e-6f, "bound ordering");
    check(r.lower_bound < 1e-4f, "identical meshes have zero sampled distance");
    check(r.upper_bound < 2.0f, "identical meshes have a tight error bound");

    /* Translate one box by 0.25 along x. The true Hausdorff distance is 0.25. */
    for (size_t i = 0; i < b.vertex_count; ++i) {
        b.vertices[i].x += 0.25f;
    }
    check(geomlib_hausdorff_bounded(&a, &b, &r) == GEOMLIB_MESH_SUCCESS, "translated rc");
    check(std::fabs(r.lower_bound - 0.25f) < 1e-3f, "translated lower bound");
    check(r.upper_bound >= 0.25f, "translated upper bound contains truth");
    check(r.upper_bound - r.lower_bound < 1.8f, "upper bound is reasonably tight");

    geomlib_trimesh_free(&a);
    geomlib_trimesh_free(&b);
    if (failures != 0) {
        std::printf("%d distance test(s) failed\n", failures);
        return 1;
    }
    std::printf("All bounded Hausdorff tests passed\n");
    return 0;
}
