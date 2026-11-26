/*
 * 2D segment intersection: crossing, endpoint touching, collinear overlap and
 * disjoint cases.
 */
#include <geomlib/geomlib_segment2d.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_seg2d_result r=geomlib_segment2d_intersect(
        geomlib_vec2_create(-1,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(0,-1),geomlib_vec2_create(0,1));
    std::printf("cross kind=%d point=(%.2f,%.2f) s=%.2f t=%.2f\n",
                (int)r.kind,r.point.x,r.point.y,r.s,r.t);
    check(r.kind==GEOMLIB_SEG2D_PROPER && std::fabs(r.point.x)<1e-5 &&
          std::fabs(r.s-0.5)<1e-5,"crossing at origin");

    r=geomlib_segment2d_intersect(
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,0),geomlib_vec2_create(2,1));
    check(r.kind==GEOMLIB_SEG2D_ENDPOINT,"endpoint touch");

    r=geomlib_segment2d_intersect(
        geomlib_vec2_create(0,0),geomlib_vec2_create(2,0),
        geomlib_vec2_create(1,0),geomlib_vec2_create(3,0));
    std::printf("collinear kind=%d overlap x=(%.2f,%.2f)\n",
                (int)r.kind,r.overlap[0].x,r.overlap[1].x);
    check(r.kind==GEOMLIB_SEG2D_COLLINEAR &&
          std::fabs(r.overlap[0].x-1.0)<1e-5 &&
          std::fabs(r.overlap[1].x-2.0)<1e-5,"collinear overlap [1,2]");

    r=geomlib_segment2d_intersect(
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(0,2),geomlib_vec2_create(1,2));
    check(r.kind==GEOMLIB_SEG2D_NONE,"parallel disjoint");
    if(failures!=0){std::printf("%d segment2d test(s) failed\n",failures);return 1;}
    std::printf("All segment2d tests passed\n");
    return 0;
}
