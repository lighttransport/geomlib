/*
 * Exact minimum enclosing circle: unit circle, two points, and an obtuse
 * triangle (whose smallest circle is the longest-edge diameter).
 */
#include <geomlib/geomlib_min_circle2d.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static bool encloses(const std::vector<geomlib_vec2> &p,
                     const geomlib_min_circle2d &c)
{
    for(auto &q:p)
        if(std::hypot((double)q.x-c.center.x,(double)q.y-c.center.y)>c.radius+1e-3)
            return false;
    return true;
}

int main(void)
{
    std::vector<geomlib_vec2> circle;
    for(int i=0;i<200;++i){
        float a=i*6.283185f/200.0f;
        circle.push_back(geomlib_vec2_create(std::cos(a),std::sin(a)));
    }
    geomlib_min_circle2d c;
    int rc=geomlib_minimum_enclosing_circle(circle.data(),circle.size(),&c);
    std::printf("mincircle circle rc=%d r=%.6f\n",rc,c.radius);
    check(rc==GEOMLIB_MESH_SUCCESS,"circle computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(c.radius-1.0)<1e-3,"unit circle radius is 1");
        check(encloses(circle,c),"all points enclosed");
    }

    geomlib_vec2 two[2]={geomlib_vec2_create(-2,0),geomlib_vec2_create(2,0)};
    rc=geomlib_minimum_enclosing_circle(two,2,&c);
    std::printf("mincircle two rc=%d r=%.4f\n",rc,c.radius);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(c.radius-2.0)<1e-4,
          "two-point circle is half the distance");

    geomlib_vec2 obt[3]={geomlib_vec2_create(0,0),geomlib_vec2_create(4,0),
                         geomlib_vec2_create(1,0.5f)};
    rc=geomlib_minimum_enclosing_circle(obt,3,&c);
    std::printf("mincircle obtuse rc=%d r=%.4f\n",rc,c.radius);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(c.radius-2.0)<1e-3,
          "obtuse triangle uses longest-edge diameter");
    if(failures!=0){std::printf("%d min circle test(s) failed\n",failures);return 1;}
    std::printf("All min circle tests passed\n");
    return 0;
}
