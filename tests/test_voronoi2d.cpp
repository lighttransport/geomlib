/*
 * Voronoi tests: cells tile the bounding box, vertices are equidistant to
 * their site and no closer to another site, and locate agrees with containment.
 */
#include <geomlib/geomlib_voronoi2d.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static bool in_cell(const geomlib_voronoi2d &v, size_t i, geomlib_vec2 q)
{
    size_t a = v.offsets[i], b = v.offsets[i+1];
    if (b < a+3) return false;
    for (size_t n = a; n < b; ++n) {
        size_t m = (n+1 < b) ? n+1 : a;
        double ax=v.cells[n].x, ay=v.cells[n].y;
        double bx=v.cells[m].x, by=v.cells[m].y;
        double cross = (bx-ax)*((double)q.y-ay) - (by-ay)*((double)q.x-ax);
        if (cross < -1e-3) return false;
    }
    return true;
}

int main(void)
{
    std::vector<geomlib_vec2> pts;
    for (int gy = 0; gy < 5; ++gy)
        for (int gx = 0; gx < 5; ++gx)
            pts.push_back(geomlib_vec2_create(
                0.15f + 0.18f*gx, 0.15f + 0.18f*gy));
    geomlib_vec2 bmin = geomlib_vec2_create(0.0f, 0.0f);
    geomlib_vec2 bmax = geomlib_vec2_create(1.0f, 1.0f);

    geomlib_voronoi2d v = geomlib_voronoi2d_create();
    int rc = geomlib_compute_voronoi2d(pts.data(), pts.size(), bmin, bmax, &v);
    std::printf("voronoi rc=%d sites=%zu vertices=%zu\n",
                rc, v.site_count, v.vertex_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "voronoi builds");
    if (rc != GEOMLIB_MESH_SUCCESS) return 1;
    check(v.site_count == pts.size(), "one cell per site");

    double area = 0.0;
    for (size_t i = 0; i < v.site_count; ++i) area += geomlib_voronoi2d_cell_area(&v, i);
    std::printf("voronoi total area=%.6f (box=1)\n", area);
    check(std::fabs(area-1.0) < 1e-3, "cells tile the bounding box");

    int bad = 0, mismatch = 0, tested = 0;
    for (size_t i = 0; i < v.site_count; ++i) {
        size_t a = v.offsets[i], b = v.offsets[i+1];
        for (size_t n = a; n < b; ++n) {
            geomlib_vec2 p = v.cells[n];
            double dself = std::hypot((double)p.x-v.sites[i].x,
                                      (double)p.y-v.sites[i].y);
            for (size_t j = 0; j < v.site_count; ++j) {
                if (j == i) continue;
                double dj = std::hypot((double)p.x-v.sites[j].x,
                                       (double)p.y-v.sites[j].y);
                if (dj < dself-2e-3) bad++;
            }
        }
    }
    for (int k = 0; k < 400; ++k) {
        geomlib_vec2 q = geomlib_vec2_create(
            ((k*37)%101)/100.0f, ((k*73)%97)/100.0f);
        size_t near = geomlib_voronoi2d_locate(&v, q);
        int contained = -1;
        for (size_t i = 0; i < v.site_count; ++i)
            if (in_cell(v, i, q)) { contained = (int)i; break; }
        if (contained >= 0) { tested++; if ((size_t)contained != near) mismatch++; }
    }
    std::printf("voronoi vertex violations=%d locate mismatches=%d/%d\n",
                bad, mismatch, tested);
    check(bad == 0, "cell vertices are nearest to their own site");
    check(mismatch == 0, "locate agrees with cell containment");

    geomlib_voronoi2d_free(&v);
    if (failures != 0) {
        std::printf("%d voronoi test(s) failed\n", failures);
        return 1;
    }
    std::printf("All Voronoi tests passed\n");
    return 0;
}
