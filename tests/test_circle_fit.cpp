/*
 * Least-squares circle fitting recovers centre and radius of noisy samples.
 */
#include <geomlib/geomlib_circle_fit.h>

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
    std::vector<geomlib_vec2> pts;
    for(int i=0;i<200;++i){
        float a=i*6.283185f/200.0f;
        pts.push_back(geomlib_vec2_create(2.0f+3.0f*std::cos(a),
                                           -1.0f+3.0f*std::sin(a)));
    }
    geomlib_circle_fit fit;
    int rc=geomlib_fit_circle(pts.data(),pts.size(),&fit);
    std::printf("fit circle rc=%d center=(%.4f,%.4f) radius=%.4f rms=%.6f\n",
                rc,fit.center.x,fit.center.y,fit.radius,fit.rms);
    check(rc==GEOMLIB_MESH_SUCCESS,"circle fit computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(fit.center.x-2.0f)<1e-3f &&
              std::fabs(fit.center.y+1.0f)<1e-3f,"centre recovered");
        check(std::fabs(fit.radius-3.0f)<1e-3f,"radius recovered");
        check(fit.rms<1e-3f,"low residual");
    }
    if(failures!=0){std::printf("%d circle fit test(s) failed\n",failures);return 1;}
    std::printf("All circle fit tests passed\n");
    return 0;
}
