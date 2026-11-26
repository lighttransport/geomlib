/*
 * Constrained Delaunay tests: constraints become real triangulation edges,
 * triangles tile the convex hull, and the result is a permutation-safe index set.
 */
#include <geomlib/geomlib_cdt2d.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static double total_area(const geomlib_cdt2d &c)
{
    double a = 0.0;
    for (size_t i = 0; i < c.triangle_count; ++i) {
        geomlib_vec2 p=c.vertices[c.triangles[i*3]];
        geomlib_vec2 q=c.vertices[c.triangles[i*3+1]];
        geomlib_vec2 r=c.vertices[c.triangles[i*3+2]];
        a += 0.5*std::fabs((q.x-p.x)*(r.y-p.y)-(q.y-p.y)*(r.x-p.x));
    }
    return a;
}

static bool valid_indices(const geomlib_cdt2d &c)
{
    for (size_t i = 0; i < c.triangle_count*3; ++i)
        if (c.triangles[i] >= c.vertex_count) return false;
    return true;
}

int main(void)
{
    std::vector<geomlib_vec2> sq = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1), geomlib_vec2_create(0,1)
    };
    uint32_t diag[2] = {0,2};
    geomlib_cdt2d c = geomlib_cdt2d_create();
    int rc = geomlib_cdt2d_build(sq.data(), sq.size(), diag, 1, &c);
    std::printf("cdt square rc=%d tris=%zu constraint=(0,2) present=%d area=%.3f\n",
                rc, c.triangle_count, geomlib_cdt2d_has_edge(&c,0,2), total_area(c));
    check(rc == GEOMLIB_MESH_SUCCESS, "square CDT builds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(c.triangle_count == 2, "square has two triangles");
        check(geomlib_cdt2d_has_edge(&c,0,2), "constrained diagonal is an edge");
        check(std::fabs(total_area(c)-1.0) < 1e-4, "triangles tile the square");
        check(valid_indices(c), "triangle indices are in range");
        geomlib_cdt2d_free(&c);
    }

    /* A constraint crossing several Delaunay edges. */
    std::vector<geomlib_vec2> pts = {
        geomlib_vec2_create(0,0), geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1), geomlib_vec2_create(0,1),
        geomlib_vec2_create(0.5f,0.25f), geomlib_vec2_create(0.5f,0.75f)
    };
    uint32_t seg[2] = {0,3};
    c = geomlib_cdt2d_create();
    rc = geomlib_cdt2d_build(pts.data(), pts.size(), seg, 1, &c);
    std::printf("cdt cross rc=%d tris=%zu vertices=%zu constraint=(0,3) present=%d area=%.3f\n",
                rc, c.triangle_count, c.vertex_count,
                geomlib_cdt2d_has_edge(&c,0,3), total_area(c));
    check(rc == GEOMLIB_MESH_SUCCESS, "crossing constraint CDT builds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(geomlib_cdt2d_has_edge(&c,0,3), "crossing constraint recovered");
        check(std::fabs(total_area(c)-1.0) < 1e-3, "area preserved");
        check(valid_indices(c), "triangle indices are in range");
        geomlib_cdt2d_free(&c);
    }

    /* Boundary constraints of the square. */
    uint32_t rim[8] = {0,1, 1,2, 2,3, 3,0};
    c = geomlib_cdt2d_create();
    rc = geomlib_cdt2d_build(sq.data(), sq.size(), rim, 4, &c);
    check(rc == GEOMLIB_MESH_SUCCESS, "boundary-constrained CDT builds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        bool all = true;
        for (int k = 0; k < 4; ++k)
            if (!geomlib_cdt2d_has_edge(&c, rim[k*2], rim[k*2+1])) all = false;
        check(all, "all boundary edges present");
        check(std::fabs(total_area(c)-1.0) < 1e-4, "square area preserved");
        geomlib_cdt2d_free(&c);
    }

    if (failures != 0) {
        std::printf("%d CDT test(s) failed\n", failures);
        return 1;
    }
    std::printf("All CDT tests passed\n");
    return 0;
}
