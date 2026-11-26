/*
 * Polygon offset tests: a unit square offsets to the expected areas for
 * outward and inward distances.
 */
#include <geomlib/geomlib_polygon_offset.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static geomlib_polygon2d square()
{
    geomlib_polygon2d p;
    geomlib_polygon2d_init(&p, 4);
    geomlib_polygon2d_add_vertex(&p, geomlib_vec2_create(0,0));
    geomlib_polygon2d_add_vertex(&p, geomlib_vec2_create(1,0));
    geomlib_polygon2d_add_vertex(&p, geomlib_vec2_create(1,1));
    geomlib_polygon2d_add_vertex(&p, geomlib_vec2_create(0,1));
    return p;
}

int main(void)
{
    geomlib_polygon2d p = square();
    geomlib_offset_result out = geomlib_offset_result_create();
    int rc = geomlib_offset_polygon(&p, 0.1f, 0, &out);
    std::printf("outward rc=%d loops=%zu\n", rc, out.count);
    check(rc == GEOMLIB_MESH_SUCCESS && out.count == 1, "outward offset succeeds");
    if (out.count == 1) {
        float area = geomlib_polygon2d_area(&out.polygons[0]);
        std::printf("  area=%.4f expected=1.44\n", area);
        check(std::fabs(fabs(area) - 1.44f) < 0.01f, "outward area is 1.44");
    }
    geomlib_offset_result_free(&out);

    out = geomlib_offset_result_create();
    rc = geomlib_offset_polygon(&p, -0.1f, 0, &out);
    std::printf("inward rc=%d loops=%zu\n", rc, out.count);
    check(rc == GEOMLIB_MESH_SUCCESS && out.count == 1, "inward offset succeeds");
    if (out.count == 1) {
        float area = geomlib_polygon2d_area(&out.polygons[0]);
        std::printf("  area=%.4f expected=0.64\n", area);
        check(std::fabs(fabs(area) - 0.64f) < 0.01f, "inward area is 0.64");
    }
    geomlib_offset_result_free(&out);
    geomlib_polygon2d_free(&p);

    if (failures != 0) {
        std::printf("%d polygon offset test(s) failed\n", failures);
        return 1;
    }
    std::printf("All polygon offset tests passed\n");
    return 0;
}
