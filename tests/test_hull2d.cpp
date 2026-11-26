/*
 * 2D convex hull tests: square with interior points, collinear inputs,
 * CCW convexity and containment of every input point.
 */
#include <geomlib/geomlib_hull2d.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static double cross(const geomlib_vec2 &o, const geomlib_vec2 &a, const geomlib_vec2 &b)
{
    return (a.x-o.x)*(b.y-o.y)-(a.y-o.y)*(b.x-o.x);
}

int main(void)
{
    std::vector<geomlib_vec2> pts = {
        geomlib_vec2_create(-1,-1), geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1), geomlib_vec2_create(-1,1),
        geomlib_vec2_create(0,0), geomlib_vec2_create(0.3f,-0.2f),
        geomlib_vec2_create(-1,0), geomlib_vec2_create(0,1)
    };
    geomlib_hull2d h = geomlib_hull2d_create();
    int rc = geomlib_convex_hull2d(pts.data(), pts.size(), &h);
    std::printf("hull square rc=%d vertices=%zu\n", rc, h.count);
    check(rc == GEOMLIB_MESH_SUCCESS, "square hull builds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(h.count == 4, "square hull has four vertices");
        bool convex = true;
        for (size_t i = 0; i < h.count; ++i)
            if (cross(h.vertices[i], h.vertices[(i+1)%h.count],
                      h.vertices[(i+2)%h.count]) < 1e-6) convex = false;
        check(convex, "hull is CCW and strictly convex");
        bool contains = true;
        for (auto &p : pts)
            for (size_t i = 0; i < h.count; ++i)
                if (cross(h.vertices[i], h.vertices[(i+1)%h.count], p) < -1e-5)
                    contains = false;
        check(contains, "all input points lie inside the hull");
        double area = 0.0;
        for (size_t i = 0; i < h.count; ++i) {
            size_t j = (i+1)%h.count;
            area += h.vertices[i].x*h.vertices[j].y - h.vertices[j].x*h.vertices[i].y;
        }
        std::printf("hull area=%.4f (expect 4)\n", 0.5*area);
        check(std::fabs(0.5*area-4.0) < 1e-4, "square hull area is 4");
        geomlib_hull2d_free(&h);
    }

    std::vector<geomlib_vec2> line = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(1,1),
        geomlib_vec2_create(2,2), geomlib_vec2_create(0.5f,0.5f)
    };
    rc = geomlib_convex_hull2d(line.data(), line.size(), &h);
    std::printf("hull collinear rc=%d vertices=%zu\n", rc, h.count);
    check(rc == GEOMLIB_MESH_SUCCESS && h.count == 2,
          "collinear hull has two endpoints");
    geomlib_hull2d_free(&h);

    if (failures != 0) {
        std::printf("%d hull test(s) failed\n", failures);
        return 1;
    }
    std::printf("All hull tests passed\n");
    return 0;
}
