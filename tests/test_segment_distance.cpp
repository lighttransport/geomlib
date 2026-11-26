/*
 * Segment distance: crossing segments give zero, parallel segments give the
 * separation, and skew segments give a positive distance with valid params.
 */
#include <geomlib/geomlib_segment_distance.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_segment_distance d;
    int rc=geomlib_segment_closest_points(
        geomlib_vec3_create(-1,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,-1,0),geomlib_vec3_create(0,1,0),&d);
    std::printf("cross rc=%d d=%.5f s=%.2f t=%.2f\n",rc,d.distance,d.s,d.t);
    check(rc==GEOMLIB_MESH_SUCCESS && d.distance<1e-5f,"crossing segments touch");

    rc=geomlib_segment_closest_points(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,2,0),geomlib_vec3_create(1,2,0),&d);
    std::printf("parallel rc=%d d=%.5f\n",rc,d.distance);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(d.distance-2.0f)<1e-4f,
          "parallel separation is 2");

    rc=geomlib_segment_closest_points(
        geomlib_vec3_create(-1,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,0,1),geomlib_vec3_create(0,0,3),&d);
    std::printf("skew rc=%d d=%.5f s=%.2f t=%.2f\n",rc,d.distance,d.s,d.t);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(d.distance-1.0f)<1e-4f &&
          d.s>=0.0f && d.s<=1.0f && d.t>=0.0f && d.t<=1.0f,
          "skew distance is 1 with valid parameters");
    if(failures!=0){std::printf("%d segment distance test(s) failed\n",failures);return 1;}
    std::printf("All segment distance tests passed\n");
    return 0;
}
