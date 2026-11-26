/*
 * Minkowski sum tests: square + square area, commutativity, and translation
 * by a single point.
 */
#include <geomlib/geomlib_minkowski2d.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static double area(const geomlib_hull2d &h)
{
    double s = 0.0;
    for (size_t i = 0; i < h.count; ++i) {
        size_t j=(i+1)%h.count;
        s += (double)h.vertices[i].x*h.vertices[j].y
           - (double)h.vertices[j].x*h.vertices[i].y;
    }
    return 0.5*std::fabs(s);
}

static bool convex(const geomlib_hull2d &h)
{
    for (size_t i = 0; i < h.count; ++i) {
        const auto &o=h.vertices[i], &a=h.vertices[(i+1)%h.count], &b=h.vertices[(i+2)%h.count];
        if ((a.x-o.x)*(b.y-o.y)-(a.y-o.y)*(b.x-o.x) < -1e-6) return false;
    }
    return true;
}

int main(void)
{
    std::vector<geomlib_vec2> sq = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1), geomlib_vec2_create(0,1)
    };
    geomlib_vec2 origin = geomlib_vec2_create(2,3);

    geomlib_hull2d h = geomlib_hull2d_create();
    int rc = geomlib_minkowski_sum2d(sq.data(), sq.size(), sq.data(), sq.size(), &h);
    std::printf("minkowski square+square rc=%d vertices=%zu area=%.4f\n", rc, h.count, area(h));
    check(rc == GEOMLIB_MESH_SUCCESS, "square sum builds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(h.count == 4 && std::fabs(area(h)-4.0) < 1e-4, "sum is a 2x2 square");
        check(convex(h), "sum is convex");
        geomlib_hull2d_free(&h);
    }

    geomlib_hull2d h2 = geomlib_hull2d_create();
    rc = geomlib_minkowski_sum2d(sq.data(), sq.size(), &origin, 1, &h2);
    std::printf("minkowski square+point rc=%d vertices=%zu area=%.4f\n", rc, h2.count, area(h2));
    check(rc == GEOMLIB_MESH_SUCCESS && h2.count == 4 &&
          std::fabs(area(h2)-1.0) < 1e-4, "point sum translates the square");
    geomlib_hull2d_free(&h2);

    geomlib_hull2d h3 = geomlib_hull2d_create();
    std::vector<geomlib_vec2> tri = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(1,0), geomlib_vec2_create(0,1)
    };
    rc = geomlib_minkowski_sum2d(tri.data(), tri.size(), sq.data(), sq.size(), &h3);
    std::printf("minkowski triangle+square rc=%d vertices=%zu area=%.4f\n", rc, h3.count, area(h3));
    check(rc == GEOMLIB_MESH_SUCCESS && h3.count >= 4 && convex(h3),
          "triangle sum is convex");
    geomlib_hull2d_free(&h3);

    if (failures != 0) {
        std::printf("%d Minkowski test(s) failed\n", failures);
        return 1;
    }
    std::printf("All Minkowski tests passed\n");
    return 0;
}
