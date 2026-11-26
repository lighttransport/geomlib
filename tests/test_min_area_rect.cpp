/*
 * Minimum-area rectangle: squares and rotated rectangles give analytic area,
 * and every input point lies inside the returned box.
 */
#include <geomlib/geomlib_min_area_rect.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static bool inside(const geomlib_min_area_rect &r, geomlib_vec2 p)
{
    double dx=(double)p.x-r.center.x,dy=(double)p.y-r.center.y;
    for(int k=0;k<2;++k){
        double proj=dx*r.axes[k].x+dy*r.axes[k].y;
        if(std::fabs(proj)>((k==0)?r.extents.x:r.extents.y)+1e-3) return false;
    }
    return true;
}

int main(void)
{
    std::vector<geomlib_vec2> sq={
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1),geomlib_vec2_create(0,1)};
    geomlib_min_area_rect r;
    int rc=geomlib_minimum_area_rectangle(sq.data(),sq.size(),&r);
    std::printf("rect square rc=%d area=%.4f extents=(%.3f,%.3f)\n",
                rc,r.area,r.extents.x,r.extents.y);
    check(rc==GEOMLIB_MESH_SUCCESS,"square rectangle computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(r.area-1.0)<1e-4,"square area is 1");
        check(std::fabs(r.center.x-0.5)<1e-3 && std::fabs(r.center.y-0.5)<1e-3,
              "square centre is correct");
        bool in=true; for(auto &p:sq) if(!inside(r,p)) in=false;
        check(in,"all points inside rectangle");
    }

    std::vector<geomlib_vec2> rot;
    for(int i=0;i<4;++i){
        float x=(i==1||i==2)?2.0f:0.0f;
        float y=(i>=2)?1.0f:0.0f;
        float a=0.5235988f;
        rot.push_back(geomlib_vec2_create(
            x*std::cos(a)-y*std::sin(a),x*std::sin(a)+y*std::cos(a)));
    }
    rc=geomlib_minimum_area_rectangle(rot.data(),rot.size(),&r);
    std::printf("rect rotated rc=%d area=%.4f expect=2\n",rc,r.area);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(r.area-2.0)<1e-3,
          "rotated rectangle area is 2");
    if(rc==GEOMLIB_MESH_SUCCESS){
        bool in=true; for(auto &p:rot) if(!inside(r,p)) in=false;
        check(in,"all rotated points inside rectangle");
    }

    std::vector<geomlib_vec2> tri={
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(0,1)};
    rc=geomlib_minimum_area_rectangle(tri.data(),tri.size(),&r);
    std::printf("rect triangle rc=%d area=%.4f\n",rc,r.area);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(r.area-1.0)<1e-3,
          "right triangle rectangle area is 1");
    if(failures!=0){std::printf("%d min rect test(s) failed\n",failures);return 1;}
    std::printf("All min rect tests passed\n");
    return 0;
}
