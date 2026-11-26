/*
 * Dynamic BVH refit test: after translating a closed mesh, refitted closest
 * point queries agree with a freshly rebuilt tree.
 */
#include <geomlib/geomlib_bvh_refit.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh mesh, moved;
    geomlib_create_sphere(&mesh, 1.0f, 16, 8);
    geomlib_trimesh_init(&moved, mesh.vertex_count, mesh.triangle_count);
    for (size_t i = 0; i < mesh.vertex_count; ++i) {
        geomlib_vec3 p = mesh.vertices[i];
        p.x += 3.0f;
        geomlib_trimesh_add_vertex(&moved, p, nullptr, nullptr);
    }
    for (size_t i = 0; i < mesh.triangle_count; ++i)
        geomlib_trimesh_add_triangle(&moved, mesh.indices[i*3],
                                     mesh.indices[i*3+1], mesh.indices[i*3+2]);

    geomlib_bvh bvh, fresh;
    geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);
    geomlib_bvh_refit(&bvh, &moved);
    geomlib_bvh_build(&fresh, &moved, GEOMLIB_BVH_QUALITY_MEDIUM);

    geomlib_vec3 q = geomlib_vec3_create(3.0f, 0.0f, 1.5f);
    geomlib_vec3 p1, p2;
    float d1, d2;
    bool h1 = geomlib_bvh_closest_point(&bvh, q, &p1, &d1, NULL);
    bool h2 = geomlib_bvh_closest_point(&fresh, q, &p2, &d2, NULL);
    std::printf("refit hit=%d d=%.4f fresh hit=%d d=%.4f\n", h1, d1, h2, d2);
    check(h1 && h2, "both trees answer closest-point queries");
    check(std::fabs(d1-d2) < 1e-3f, "refit matches a rebuilt tree");
    check(p1.x > 2.0f, "refit bounds follow the moved mesh");

    geomlib_bvh_free(&bvh);
    geomlib_bvh_free(&fresh);
    geomlib_trimesh_free(&moved);
    geomlib_trimesh_free(&mesh);
    if (failures != 0) {
        std::printf("%d BVH refit test(s) failed\n", failures);
        return 1;
    }
    std::printf("All BVH refit tests passed\n");
    return 0;
}
