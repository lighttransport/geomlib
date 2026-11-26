/*
 * Medial axis tests: sampled Voronoi centres lie inside the polygon and are
 * equidistant to at least two boundary edges (centres of maximal circles).
 */
#include <geomlib/geomlib_medial2d.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static bool inside(geomlib_vec2 q, const geomlib_vec2 *p, size_t n)
{
    bool in = false;
    for (size_t i = 0; i < n; ++i) {
        size_t j = (i+1)%n;
        if (((p[i].y > q.y) != (p[j].y > q.y)) &&
            ((double)q.x < ((double)p[j].x-p[i].x)*((double)q.y-p[i].y)/
                           ((double)p[j].y-p[i].y)+p[i].x))
            in = !in;
    }
    return in;
}

static double seg_dist(geomlib_vec2 q, geomlib_vec2 a, geomlib_vec2 b)
{
    double vx=(double)b.x-a.x, vy=(double)b.y-a.y;
    double wx=(double)q.x-a.x, wy=(double)q.y-a.y;
    double l2=vx*vx+vy*vy, t = l2>0 ? (wx*vx+wy*vy)/l2 : 0.0;
    if (t<0) t=0; if (t>1) t=1;
    return std::hypot((double)q.x-(a.x+t*(b.x-a.x)),
                      (double)q.y-(a.y+t*(b.y-a.y)));
}

int main(void)
{
    geomlib_vec2 square[4] = {
        geomlib_vec2_create(-1,-1), geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1), geomlib_vec2_create(-1,1)
    };
    geomlib_medial2d m = geomlib_medial2d_create();
    int rc = geomlib_medial_axis2d(square, 4, 0.06f, &m);
    std::printf("medial square rc=%d segments=%zu\n", rc, m.segment_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "square medial axis builds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(m.segment_count > 0, "medial axis has segments");
        int outside = 0, not_two = 0;
        double maxr = 0.0;
        for (size_t i = 0; i < m.segment_count*2; ++i) {
            geomlib_vec2 q = m.segments[i];
            if (!inside(q, square, 4)) outside++;
            double d[4];
            for (int k = 0; k < 4; ++k)
                d[k] = seg_dist(q, square[k], square[(k+1)%4]);
            std::sort(d, d+4);
            if (d[1] - d[0] > 0.06) not_two++;
            if (d[0] > maxr) maxr = d[0];
        }
        std::printf("outside=%d not_two=%d max_inradius=%.3f\n",
                    outside, not_two, maxr);
        check(outside == 0, "medial endpoints inside polygon");
        check(not_two == 0, "medial endpoints equidistant to two edges");
        check(maxr > 0.5, "medial axis reaches the incenter");
        geomlib_medial2d_free(&m);
    }

    geomlib_vec2 rect[4] = {
        geomlib_vec2_create(-2,-1), geomlib_vec2_create(2,-1),
        geomlib_vec2_create(2,1), geomlib_vec2_create(-2,1)
    };
    m = geomlib_medial2d_create();
    rc = geomlib_medial_axis2d(rect, 4, 0.08f, &m);
    std::printf("medial rect rc=%d segments=%zu\n", rc, m.segment_count);
    check(rc == GEOMLIB_MESH_SUCCESS && m.segment_count > 0,
          "rectangle medial axis builds");
    geomlib_medial2d_free(&m);

    if (failures != 0) {
        std::printf("%d medial axis test(s) failed\n", failures);
        return 1;
    }
    std::printf("All medial axis tests passed\n");
    return 0;
}
