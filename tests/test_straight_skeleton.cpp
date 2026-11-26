/*
 * Straight skeleton tests: a square's skeleton should trace its four corners
 * toward the centre, producing nodes and edges with finite coordinates.
 */
#include <geomlib/geomlib_straight_skeleton.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_polygon2d square;
    geomlib_polygon2d_init(&square, 4);
    geomlib_vec2 p[4] = {{0,0},{1,0},{1,1},{0,1}};
    for (int i = 0; i < 4; ++i) geomlib_polygon2d_add_vertex(&square, p[i]);

    geomlib_straight_skeleton s = geomlib_straight_skeleton_create();
    int rc = geomlib_compute_straight_skeleton(&square, &s);
    std::printf("straight skeleton rc=%d nodes=%zu edges=%zu\n",
                rc, s.point_count, s.edge_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "straight skeleton succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(s.point_count >= 4 && s.edge_count > 0, "skeleton has nodes and edges");
        bool valid = true;
        for (size_t i = 0; i < s.point_count; ++i) {
            if (!std::isfinite(s.points[i].x) || !std::isfinite(s.points[i].y) ||
                !std::isfinite(s.times[i])) valid = false;
            if (s.points[i].x < -1.0f || s.points[i].x > 2.0f ||
                s.points[i].y < -1.0f || s.points[i].y > 2.0f) valid = false;
        }
        for (size_t i = 0; i < s.edge_count*2; ++i)
            if (s.edges[i] >= s.point_count) valid = false;
        check(valid, "skeleton nodes and edges are valid and in bounds");
        /* At least one traced point should approach the centre. */
        bool near_center = false;
        for (size_t i = 0; i < s.point_count; ++i)
            if (std::fabs(s.points[i].x-0.5f) < 0.25f &&
                std::fabs(s.points[i].y-0.5f) < 0.25f) near_center = true;
        check(near_center, "offset wavefront approaches the square centre");
        geomlib_straight_skeleton_free(&s);
    }
    geomlib_polygon2d_free(&square);
    if (failures != 0) {
        std::printf("%d straight skeleton test(s) failed\n", failures);
        return 1;
    }
    std::printf("All straight skeleton tests passed\n");
    return 0;
}
