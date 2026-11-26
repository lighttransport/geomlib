/*
 * Mesh intersection tests: two overlapping boxes produce a closed intersection
 * loop, disjoint boxes produce none.
 */
#include <geomlib/geomlib_intersection.h>
#include <geomlib/geomlib_bvh.h>

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
    geomlib_trimesh a = {}, b = {};
    check(geomlib_create_box(&a, 1.0f,1.0f,1.0f) == GEOMLIB_MESH_SUCCESS, "box A");
    check(geomlib_create_box(&b, 1.0f,1.0f,1.0f) == GEOMLIB_MESH_SUCCESS, "box B");
    for (size_t i = 0; i < b.vertex_count; ++i) {
        b.vertices[i].x += 0.5f;
        b.vertices[i].y += 0.3f;
        b.vertices[i].z += 0.2f;
    }

    geomlib_intersection r = geomlib_mesh_mesh_intersection(&a, &b);
    std::printf("intersection segments=%zu\n", r.segment_count);
    check(r.segment_count > 0, "overlapping boxes intersect");

    /* Every endpoint must lie in the overlap region [-0.5,0.5] x [-0.5,0.5]?
     * The overlap is x in [0,0.5], y in [-0.5,0.5], z in [-0.5,0.5]. */
    int bad = 0;
    for (size_t i = 0; i < r.segment_count*2; ++i) {
        geomlib_vec3 p = r.segments[i];
        if (p.x < -0.01f || p.x > 0.51f ||
            p.y < -0.51f || p.y > 0.51f ||
            p.z < -0.51f || p.z > 0.51f) bad++;
    }
    std::printf("endpoints outside overlap=%d\n", bad);
    check(bad == 0, "intersection points lie in the overlap region");

    geomlib_bvh bvha, bvhb;
    geomlib_bvh_build(&bvha, &a, GEOMLIB_BVH_QUALITY_MEDIUM);
    geomlib_bvh_build(&bvhb, &b, GEOMLIB_BVH_QUALITY_MEDIUM);
    int offsurface = 0;
    for (size_t i = 0; i < r.segment_count*2; ++i) {
        float da=0.0f, db=0.0f;
        geomlib_vec3 cp;
        geomlib_bvh_closest_point(&bvha, r.segments[i], &cp, &da, nullptr);
        geomlib_bvh_closest_point(&bvhb, r.segments[i], &cp, &db, nullptr);
        if (da > 1e-3f || db > 1e-3f) offsurface++;
    }
    std::printf("endpoints off both surfaces=%d\n", offsurface);
    check(offsurface == 0, "intersection points lie on both surfaces");
    geomlib_bvh_free(&bvha);
    geomlib_bvh_free(&bvhb);
    geomlib_intersection_free(&r);

    /* Disjoint boxes. */
    for (size_t i = 0; i < b.vertex_count; ++i) b.vertices[i].x += 2.0f;
    r = geomlib_mesh_mesh_intersection(&a, &b);
    std::printf("disjoint segments=%zu\n", r.segment_count);
    check(r.segment_count == 0, "disjoint boxes do not intersect");
    geomlib_intersection_free(&r);

    geomlib_trimesh_free(&a);
    geomlib_trimesh_free(&b);
    if (failures != 0) {
        std::printf("%d intersection test(s) failed\n", failures);
        return 1;
    }
    std::printf("All intersection tests passed\n");
    return 0;
}
