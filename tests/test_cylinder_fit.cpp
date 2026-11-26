/*
 * Cylinder fitting: recover axis, centre line and radius from sampled points.
 */
#include <geomlib/geomlib_cylinder_fit.h>

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
    for(int i=0;i<300;++i){
        float t=(float)((i%30)-15)*0.2f;
        float a=(float)(i/30)*6.283185f/10.0f;
        pts.push_back(geomlib_vec3_create(2.0f+1.5f*std::cos(a),
                                          -1.0f+1.5f*std::sin(a),t));
    }
    geomlib_cylinder_fit fit;
    int rc=geomlib_fit_cylinder(pts.data(),pts.size(),
                                geomlib_vec3_create(0,0,1),&fit);
    std::printf("fit cylinder rc=%d axis=(%.3f,%.3f,%.3f) point=(%.3f,%.3f) r=%.3f rms=%.5f\n",
                rc,fit.axis.x,fit.axis.y,fit.axis.z,fit.point.x,fit.point.y,
                fit.radius,fit.rms);
    check(rc==GEOMLIB_MESH_SUCCESS,"cylinder fit computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(fit.axis.z)>0.99f,"axis is vertical");
        check(std::fabs(fit.point.x-2.0f)<0.05f &&
              std::fabs(fit.point.y+1.0f)<0.05f,"axis centre recovered");
        check(std::fabs(fit.radius-1.5f)<0.05f,"radius recovered");
        check(fit.rms<0.05f,"low residual");
    }
    if(failures!=0){std::printf("%d cylinder fit test(s) failed\n",failures);return 1;}
    std::printf("All cylinder fit tests passed\n");
    return 0;
}
