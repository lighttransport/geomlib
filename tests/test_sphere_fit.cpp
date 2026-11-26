/*
 * Least-squares sphere fitting recovers centre and radius of sampled points.
 */
#include <geomlib/geomlib_sphere_fit.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    std::vector<geomlib_vec3> pts;
    for(int i=0;i<400;++i){
        float u=i*0.618033f*6.283185f;
        float z=1.0f-2.0f*(float)(i+0.5)/400.0f;
        float r=std::sqrt(std::fmax(0.0f,1.0f-z*z));
        pts.push_back(geomlib_vec3_create(1.0f+2.0f*r*std::cos(u),
                                          -2.0f+2.0f*r*std::sin(u),
                                          3.0f+2.0f*z));
    }
    geomlib_sphere_fit fit;
    int rc=geomlib_fit_sphere(pts.data(),pts.size(),&fit);
    std::printf("fit sphere rc=%d center=(%.4f,%.4f,%.4f) radius=%.4f rms=%.6f\n",
                rc,fit.center.x,fit.center.y,fit.center.z,fit.radius,fit.rms);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere fit computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(fit.center.x-1.0f)<1e-2f &&
              std::fabs(fit.center.y+2.0f)<1e-2f &&
              std::fabs(fit.center.z-3.0f)<1e-2f,"centre recovered");
        check(std::fabs(fit.radius-2.0f)<1e-2f,"radius recovered");
        check(fit.rms<1e-2f,"low residual");
    }
    if(failures!=0){std::printf("%d sphere fit test(s) failed\n",failures);return 1;}
    std::printf("All sphere fit tests passed\n");
    return 0;
}
