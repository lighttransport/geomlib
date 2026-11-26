/*
 * Ray/plane intersection: analytic hits, misses and segment bounds.
 */
#include <geomlib/geomlib_ray_plane.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_vec3 p0=geomlib_vec3_create(0,0,2);
    geomlib_vec3 n=geomlib_vec3_create(0,0,1);
    float t; geomlib_vec3 hit;
    int ok=geomlib_ray_plane_intersection(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(0,0,1),p0,n,&t,&hit);
    std::printf("ray hit ok=%d t=%.3f z=%.3f\n",ok,t,hit.z);
    check(ok && std::fabs(t-2.0f)<1e-4f && std::fabs(hit.z-2.0f)<1e-4f,
          "ray hits plane z=2 at t=2");

    ok=geomlib_ray_plane_intersection(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(0,0,-1),p0,n,&t,&hit);
    check(!ok,"ray pointing away misses");
    ok=geomlib_ray_plane_intersection(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),p0,n,&t,&hit);
    check(!ok,"parallel ray misses");

    ok=geomlib_segment_plane_intersection(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(0,0,4),p0,n,&t,&hit);
    check(ok && std::fabs(t-0.5f)<1e-4f,"segment crosses at t=0.5");
    ok=geomlib_segment_plane_intersection(
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(0,0,1),p0,n,&t,&hit);
    check(!ok,"short segment does not reach plane");
    if(failures!=0){std::printf("%d ray-plane test(s) failed\n",failures);return 1;}
    std::printf("All ray-plane tests passed\n");
    return 0;
}
