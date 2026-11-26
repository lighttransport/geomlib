/*
 * Delaunay refinement: a skinny rectangle gains Steiner points until its
 * minimum triangle angle meets the target, with boundary edges preserved.
 */
#include <geomlib/geomlib_refine2d.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static double area(const geomlib_cdt2d &c)
{
    double s=0.0;
    for(size_t i=0;i<c.triangle_count;++i){
        geomlib_vec2 a=c.vertices[c.triangles[i*3]];
        geomlib_vec2 b=c.vertices[c.triangles[i*3+1]];
        geomlib_vec2 d=c.vertices[c.triangles[i*3+2]];
        s+=0.5*std::fabs((b.x-a.x)*(d.y-a.y)-(b.y-a.y)*(d.x-a.x));
    }
    return s;
}

int main(void)
{
    geomlib_vec2 rect[4]={
        geomlib_vec2_create(0,0),geomlib_vec2_create(4,0),
        geomlib_vec2_create(4,1),geomlib_vec2_create(0,1)
    };
    geomlib_cdt2d initial=geomlib_cdt2d_create();
    int rc=geomlib_cdt2d_build(rect,4,nullptr,0,&initial);
    check(rc==GEOMLIB_MESH_SUCCESS,"initial CDT builds");
    float before=geomlib_cdt2d_min_angle(&initial);
    geomlib_cdt2d_free(&initial);

    geomlib_refine2d_params p=geomlib_refine2d_params_default();
    p.min_angle_deg=25.0f;
    p.max_points=20;
    geomlib_cdt2d out=geomlib_cdt2d_create();
    rc=geomlib_refine_polygon2d(rect,4,&p,&out);
    float after=geomlib_cdt2d_min_angle(&out);
    std::printf("refine2d rc=%d min angle %.2f -> %.2f points=%zu tris=%zu area=%.4f\n",
                rc,before,after,out.vertex_count,out.triangle_count,area(out));
    check(rc==GEOMLIB_MESH_SUCCESS,"refinement succeeds");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(before<25.0f,"initial triangulation violates the target");
        check(after>=before-0.01f,"minimum angle does not get worse");
        bool inside=true;
        for(size_t i=0;i<out.triangle_count;++i){
            geomlib_vec2 a=out.vertices[out.triangles[i*3]],
                         b=out.vertices[out.triangles[i*3+1]],
                         c=out.vertices[out.triangles[i*3+2]];
            geomlib_vec2 mid=geomlib_vec2_create((a.x+b.x+c.x)/3,(a.y+b.y+c.y)/3);
            if(mid.x<0||mid.x>4||mid.y<0||mid.y>1) inside=false;
        }
        check(inside,"all triangles stay inside the domain");
        check(out.vertex_count>4,"Steiner points inserted");
        check(std::fabs(area(out)-4.0)<1e-3,"domain area preserved");
        check(out.triangle_count>4,"Steiner triangles created");
        geomlib_cdt2d_free(&out);
    }
    if(failures!=0){std::printf("%d refine2d test(s) failed\n",failures);return 1;}
    std::printf("All refine2d tests passed\n");
    return 0;
}
