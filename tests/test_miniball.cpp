/*
 * Minimum enclosing sphere: on a sphere and box point sets, the exact
 * smallest enclosing sphere is recovered and is no larger than Ritter's.
 */
#include <geomlib/geomlib_miniball.h>
#include <geomlib/geomlib_bounding_sphere.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static bool encloses(const std::vector<geomlib_vec3> &p,
                     const geomlib_bounding_sphere &s)
{
    for(auto &q:p)
        if(geomlib_vec3_length(geomlib_vec3_sub(q,s.center))>s.radius+1e-3)
            return false;
    return true;
}

int main(void)
{
    std::vector<geomlib_vec3> sphere;
    for(int i=0;i<200;++i){
        float u=i*0.618033f*6.283185f;
        float z=1.0f-2.0f*(float)(i+0.5)/200.0f;
        float r=std::sqrt(std::fmax(0.0f,1.0f-z*z));
        sphere.push_back(geomlib_vec3_create(r*std::cos(u),r*std::sin(u),z));
    }
    geomlib_bounding_sphere exact,ritter;
    int rc=geomlib_minimum_enclosing_sphere(sphere.data(),sphere.size(),&exact);
    geomlib_bounding_sphere_points(sphere.data(),sphere.size(),&ritter);
    std::printf("miniball sphere rc=%d radius=%.5f ritter=%.5f\n",rc,exact.radius,ritter.radius);
    check(rc==GEOMLIB_MESH_SUCCESS,"miniball computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(encloses(sphere,exact),"all points enclosed");
        check(exact.radius>=0.99f && exact.radius<1.02f,"radius close to 1");
        check(exact.radius<=ritter.radius+1e-3f,"exact bound is no larger");
    }

    std::vector<geomlib_vec3> box;
    for(int z=-1;z<=1;z+=2)
        for(int y=-1;y<=1;y+=2)
            for(int x=-1;x<=1;x+=2)
                box.push_back(geomlib_vec3_create((float)x,(float)y,(float)z));
    rc=geomlib_minimum_enclosing_sphere(box.data(),box.size(),&exact);
    std::printf("miniball box rc=%d radius=%.5f expect=1.73205\n",rc,exact.radius);
    check(rc==GEOMLIB_MESH_SUCCESS && encloses(box,exact) &&
          std::fabs(exact.radius-1.73205f)<1e-3f,"box miniball is exact");
    if(failures!=0){std::printf("%d miniball test(s) failed\n",failures);return 1;}
    std::printf("All miniball tests passed\n");
    return 0;
}
