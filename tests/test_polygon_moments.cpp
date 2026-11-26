/*
 * Polygon moments: square and rectangle give analytic area, centroid and
 * second moments; reversing winding leaves the results unchanged.
 */
#include <geomlib/geomlib_polygon_moments.h>

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
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1),geomlib_vec2_create(0,1)
    };
    geomlib_polygon_moments m;
    int rc=geomlib_compute_polygon_moments(sq,4,&m);
    std::printf("moments square rc=%d area=%.4f c=(%.3f,%.3f) I=(%.5f,%.5f,%.5f)\n",
                rc,m.area,m.centroid.x,m.centroid.y,m.ixx,m.ixy,m.iyy);
    check(rc==GEOMLIB_MESH_SUCCESS,"square moments compute");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(m.area-1.0)<1e-9,"area is 1");
        check(std::fabs(m.centroid.x-0.5)<1e-5 && std::fabs(m.centroid.y-0.5)<1e-5,
              "centroid is (0.5,0.5)");
        check(std::fabs(m.ixx-1.0/12.0)<1e-4 && std::fabs(m.iyy-1.0/12.0)<1e-4 &&
              std::fabs(m.ixy)<1e-4,"square inertia is 1/12");
    }

    geomlib_vec2 rect[4]={
        geomlib_vec2_create(-1,-0.5),geomlib_vec2_create(1,-0.5),
        geomlib_vec2_create(1,0.5),geomlib_vec2_create(-1,0.5)
    };
    rc=geomlib_compute_polygon_moments(rect,4,&m);
    std::printf("moments rect rc=%d area=%.4f I=(%.5f,%.5f)\n",
                rc,m.area,m.ixx,m.iyy);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(m.area-2.0)<1e-9 &&
          std::fabs(m.ixx-1.0/6.0)<1e-4 && std::fabs(m.iyy-2.0/3.0)<1e-4,
          "rectangle moments are analytic");

    geomlib_vec2 rev[4]={sq[3],sq[2],sq[1],sq[0]};
    rc=geomlib_compute_polygon_moments(rev,4,&m);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(m.area-1.0)<1e-9 &&
          std::fabs(m.centroid.x-0.5)<1e-5 && std::fabs(m.ixx-1.0/12.0)<1e-4,
          "winding reversal leaves moments unchanged");
    if(failures!=0){std::printf("%d polygon moment test(s) failed\n",failures);return 1;}
    std::printf("All polygon moment tests passed\n");
    return 0;
}
