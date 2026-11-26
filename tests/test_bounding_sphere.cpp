/*
 * Ritter bounding sphere: encloses all input points with a radius close to
 * analytic for a sphere and a box.
 */
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
        if(geomlib_vec3_length(geomlib_vec3_sub(q,s.center))>s.radius+1e-4)
            return false;
    return true;
}

int main(void)
{
    std::vector<geomlib_vec3> sphere;
    for(int i=0;i<300;++i){
        float u=i*0.618033f*6.283185f;
        float z=1.0f-2.0f*(float)(i+0.5)/300.0f;
        float r=std::sqrt(std::fmax(0.0f,1.0f-z*z));
        sphere.push_back(geomlib_vec3_create(r*std::cos(u),r*std::sin(u),z));
    }
    geomlib_bounding_sphere s;
    int rc=geomlib_bounding_sphere_points(sphere.data(),sphere.size(),&s);
    std::printf("sphere radius=%.4f center=(%.3f,%.3f,%.3f)\n",
                s.radius,s.center.x,s.center.y,s.center.z);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere bound computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(encloses(sphere,s),"all sphere points enclosed");
        check(s.radius>=1.0f && s.radius<1.2f,"radius close to analytic 1");
    }

    std::vector<geomlib_vec3> box;
    for(int z=-1;z<=1;z+=2)
        for(int y=-1;y<=1;y+=2)
            for(int x=-1;x<=1;x+=2)
                box.push_back(geomlib_vec3_create((float)x,(float)y,(float)z));
    rc=geomlib_bounding_sphere_points(box.data(),box.size(),&s);
    std::printf("box radius=%.4f expect=1.7321\n",s.radius);
    check(rc==GEOMLIB_MESH_SUCCESS && encloses(box,s) &&
          std::fabs(s.radius-1.73205f)<0.05f,"box bound is tight");
    if(failures!=0){std::printf("%d bounding sphere test(s) failed\n",failures);return 1;}
    std::printf("All bounding sphere tests passed\n");
    return 0;
}
