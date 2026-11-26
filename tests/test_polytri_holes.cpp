/*
 * Polygon triangulation with holes: total triangle area equals domain area
 * and no triangle covers a hole.
 */
#include <geomlib/geomlib_polytri_holes.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static bool inside_rect(geomlib_vec2 q, float lo, float hi)
{
    return q.x > lo && q.x < hi && q.y > lo && q.y < hi;
}

int main(void)
{
    std::vector<geomlib_vec2> outer = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(4,0),
        geomlib_vec2_create(4,4), geomlib_vec2_create(0,4)
    };
    std::vector<geomlib_vec2> hole = {
        geomlib_vec2_create(1,1), geomlib_vec2_create(3,1),
        geomlib_vec2_create(3,3), geomlib_vec2_create(1,3)
    };
    const geomlib_vec2 *holes[1] = {hole.data()};
    size_t sizes[1] = {hole.size()};

    geomlib_polytri_result r = geomlib_polytri_result_create();
    int rc = geomlib_triangulate_with_holes(outer.data(), outer.size(),
                                            holes, sizes, 1, &r);
    double area = 0.0;
    int inhole = 0;
    if (rc == GEOMLIB_MESH_SUCCESS)
        for (size_t i = 0; i < r.triangle_count; ++i) {
            geomlib_vec2 a=r.vertices[r.triangles[i*3]];
            geomlib_vec2 b=r.vertices[r.triangles[i*3+1]];
            geomlib_vec2 c=r.vertices[r.triangles[i*3+2]];
            area += 0.5*std::fabs((b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x));
            geomlib_vec2 mid = geomlib_vec2_create((a.x+b.x+c.x)/3,
                                                   (a.y+b.y+c.y)/3);
            if (inside_rect(mid,1.0f,3.0f)) inhole++;
        }
    std::printf("polytri holes rc=%d tris=%zu vertices=%zu area=%.4f expect=12\n",
                rc, r.triangle_count, r.vertex_count, area);
    check(rc == GEOMLIB_MESH_SUCCESS, "hole triangulation builds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(r.triangle_count > 0, "triangles produced");
        check(std::fabs(area-12.0) < 0.05, "area excludes the hole");
        check(inhole == 0, "no triangle covers the hole");
        geomlib_polytri_result_free(&r);
    }

    /* Two disjoint holes in an 8x4 box. */
    std::vector<geomlib_vec2> outer2 = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(8,0),
        geomlib_vec2_create(8,4), geomlib_vec2_create(0,4)
    };
    std::vector<geomlib_vec2> h1 = {
        geomlib_vec2_create(1,1), geomlib_vec2_create(2,1),
        geomlib_vec2_create(2,2), geomlib_vec2_create(1,2)
    };
    std::vector<geomlib_vec2> h2 = {
        geomlib_vec2_create(5,2), geomlib_vec2_create(6,2),
        geomlib_vec2_create(6,3), geomlib_vec2_create(5,3)
    };
    const geomlib_vec2 *holes2[2] = {h1.data(), h2.data()};
    size_t sizes2[2] = {h1.size(), h2.size()};
    r = geomlib_polytri_result_create();
    rc = geomlib_triangulate_with_holes(outer2.data(), outer2.size(),
                                        holes2, sizes2, 2, &r);
    area = 0.0;
    if (rc == GEOMLIB_MESH_SUCCESS)
        for (size_t i = 0; i < r.triangle_count; ++i) {
            geomlib_vec2 a=r.vertices[r.triangles[i*3]];
            geomlib_vec2 b=r.vertices[r.triangles[i*3+1]];
            geomlib_vec2 c=r.vertices[r.triangles[i*3+2]];
            area += 0.5*std::fabs((b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x));
        }
    std::printf("polytri two holes rc=%d tris=%zu area=%.4f expect=30\n",
                rc, r.triangle_count, area);
    check(rc == GEOMLIB_MESH_SUCCESS && std::fabs(area-30.0) < 0.1,
          "two-hole area is 30");
    geomlib_polytri_result_free(&r);

    if (failures != 0) {
        std::printf("%d hole triangulation test(s) failed\n", failures);
        return 1;
    }
    std::printf("All hole triangulation tests passed\n");
    return 0;
}
