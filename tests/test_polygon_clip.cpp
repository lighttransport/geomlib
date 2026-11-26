/*
 * Convex clipping: a square clipped by a half-plane and by another square
 * gives the expected areas.
 */
#include <geomlib/geomlib_polygon_clip.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static double area(const geomlib_polygon2d &p)
{
    double s=0.0;
    for(size_t i=0;i<p.vertex_count;++i){
        size_t j=(i+1)%p.vertex_count;
        s+=(double)p.vertices[i].x*p.vertices[j].y
         -(double)p.vertices[j].x*p.vertices[i].y;
    }
    return 0.5*std::fabs(s);
}

int main(void)
{
    geomlib_vec2 sq[4]={
        geomlib_vec2_create(-1,-1),geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1),geomlib_vec2_create(-1,1)
    };
    geomlib_polygon2d out={};
    int rc=geomlib_clip_polygon_convex(sq,4,sq,4,&out);
    std::printf("clip self rc=%d vertices=%zu area=%.4f\n",rc,out.vertex_count,area(out));
    check(rc==GEOMLIB_MESH_SUCCESS && out.vertex_count==4 &&
          std::fabs(area(out)-4.0)<1e-4,"self clip reproduces the square");
    geomlib_polygon2d_free(&out);

    /* Clip with a half-plane x <= 0 represented as a large box. */
    geomlib_vec2 half[4]={
        geomlib_vec2_create(-10,-10),geomlib_vec2_create(0,-10),
        geomlib_vec2_create(0,10),geomlib_vec2_create(-10,10)
    };
    rc=geomlib_clip_polygon_convex(sq,4,half,4,&out);
    std::printf("clip half rc=%d vertices=%zu area=%.4f\n",rc,out.vertex_count,area(out));
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(area(out)-2.0)<1e-4,
          "half-plane clip gives area 2");
    geomlib_polygon2d_free(&out);

    /* Disjoint clip box. */
    geomlib_vec2 farbox[4]={
        geomlib_vec2_create(5,5),geomlib_vec2_create(6,5),
        geomlib_vec2_create(6,6),geomlib_vec2_create(5,6)
    };
    rc=geomlib_clip_polygon_convex(sq,4,farbox,4,&out);
    check(rc==GEOMLIB_MESH_SUCCESS && out.vertex_count==0,
          "disjoint clip gives empty polygon");
    geomlib_polygon2d_free(&out);
    if(failures!=0){std::printf("%d clip test(s) failed\n",failures);return 1;}
    std::printf("All clip tests passed\n");
    return 0;
}
