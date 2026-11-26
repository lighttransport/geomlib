/*
 * 2D Delaunay tests: square gives two triangles; random points satisfy the
 * empty circumcircle property.
 */
#include <geomlib/geomlib_delaunay2d.h>
#include <geomlib/geomlib_predicates.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static bool empty_circle(const geomlib_vec2 *pts,
                         const geomlib_delaunay2d_result &r, size_t t, size_t p)
{
    geomlib_vec2 a=pts[r.triangles[t*3]], b=pts[r.triangles[t*3+1]], c=pts[r.triangles[t*3+2]];
    if (p==r.triangles[t*3] || p==r.triangles[t*3+1] || p==r.triangles[t*3+2]) return true;
    int orient = geomlib_pred_orient2d(a.x,a.y,b.x,b.y,c.x,c.y);
    if (orient > 0)
        return geomlib_pred_incircle(a.x,a.y,b.x,b.y,c.x,c.y,pts[p].x,pts[p].y) <= 0;
    return geomlib_pred_incircle(a.x,a.y,c.x,c.y,b.x,b.y,pts[p].x,pts[p].y) <= 0;
}

int main(void)
{
    geomlib_vec2 square[4] = {{0,0},{1,0},{1,1},{0,1}};
    geomlib_delaunay2d_result r = geomlib_delaunay2d_result_create();
    int rc = geomlib_delaunay2d(square, 4, &r);
    std::printf("square rc=%d triangles=%zu\n", rc, r.triangle_count);
    check(rc == GEOMLIB_MESH_SUCCESS && r.triangle_count == 2,
          "square triangulates into two triangles");
    geomlib_delaunay2d_result_free(&r);

    std::vector<geomlib_vec2> pts;
    for (int i = 0; i < 40; ++i) {
        float a = i*0.618033f*6.283185f;
        float rr = 0.3f + 0.7f*((i*11)%17)/16.0f;
        pts.push_back(geomlib_vec2_create(rr*cosf(a), rr*sinf(a)));
    }
    r = geomlib_delaunay2d_result_create();
    rc = geomlib_delaunay2d(pts.data(), pts.size(), &r);
    std::printf("random rc=%d triangles=%zu\n", rc, r.triangle_count);
    check(rc == GEOMLIB_MESH_SUCCESS && r.triangle_count > 30,
          "random points produce a triangulation");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        bool empty = true;
        for (size_t t = 0; t < r.triangle_count && empty; ++t)
            for (size_t p = 0; p < pts.size(); ++p)
                if (!empty_circle(pts.data(), r, t, p)) { empty = false; break; }
        check(empty, "empty circumcircle property holds");
    }
    geomlib_delaunay2d_result_free(&r);
    if (failures != 0) {
        std::printf("%d 2D Delaunay test(s) failed\n", failures);
        return 1;
    }
    std::printf("All 2D Delaunay tests passed\n");
    return 0;
}
