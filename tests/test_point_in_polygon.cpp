/*
 * Point-in-polygon: containment, winding number and signed distance for
 * CCW and CW squares.
 */
#include <geomlib/geomlib_point_in_polygon.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_vec2 sq[4]={
        geomlib_vec2_create(-1,-1),geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1),geomlib_vec2_create(-1,1)
    };
    geomlib_vec2 q=geomlib_vec2_create(0.2f,0.3f);
    int inside=geomlib_point_in_polygon(sq,4,q);
    int wn=geomlib_polygon_winding(sq,4,q);
    float sd=geomlib_polygon_signed_distance(sq,4,q);
    std::printf("point inside=%d winding=%d sd=%.3f\n",inside,wn,sd);
    check(inside==1 && wn==1 && sd<0.0f,"interior point has negative distance");

    geomlib_vec2 out={2.0f,2.0f};
    check(geomlib_point_in_polygon(sq,4,out)==0 &&
          geomlib_polygon_winding(sq,4,out)==0 &&
          geomlib_polygon_signed_distance(sq,4,out)>0.0f,
          "exterior point has positive distance");
    check(std::fabs(geomlib_polygon_signed_distance(sq,4,out)-std::sqrt(2.0))<0.01,
          "exterior distance is analytic");

    geomlib_vec2 rev[4]={sq[3],sq[2],sq[1],sq[0]};
    check(geomlib_point_in_polygon(rev,4,q)==1 &&
          geomlib_polygon_winding(rev,4,q)==-1 &&
          geomlib_polygon_signed_distance(rev,4,q)<0.0f,
          "CW polygon still contains the point");
    if(failures!=0){std::printf("%d point-in-polygon test(s) failed\n",failures);return 1;}
    std::printf("All point-in-polygon tests passed\n");
    return 0;
}
