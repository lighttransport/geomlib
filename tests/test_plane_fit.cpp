/*
 * Least-squares plane fitting: a planar point set recovers the plane normal
 * and low residual; a tilted plane recovers its normal.
 */
#include <geomlib/geomlib_plane_fit.h>

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
    for(int i=0;i<100;++i)
        pts.push_back(geomlib_vec3_create((float)((i*37)%101)/10.0f,
                                          (float)((i*53)%97)/10.0f,1.5f));
    geomlib_plane_fit fit;
    int rc=geomlib_fit_plane(pts.data(),pts.size(),&fit);
    std::printf("fit plane rc=%d normal=(%.3f,%.3f,%.3f) rms=%.6f\n",
                rc,fit.normal.x,fit.normal.y,fit.normal.z,fit.rms);
    check(rc==GEOMLIB_MESH_SUCCESS,"plane fit computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(fit.normal.z)>0.99f,"normal is vertical");
        check(std::fabs(fit.point.z-1.5f)<1e-3f,"point lies on plane z=1.5");
        check(fit.rms<1e-4f,"low residual");
    }

    /* Tilt the points: z = 0.5 x + 1. */
    for(size_t i=0;i<pts.size();++i)
        pts[i].z=0.5f*pts[i].x+1.0f;
    rc=geomlib_fit_plane(pts.data(),pts.size(),&fit);
    std::printf("fit tilted rc=%d normal=(%.3f,%.3f,%.3f) rms=%.6f\n",
                rc,fit.normal.x,fit.normal.y,fit.normal.z,fit.rms);
    check(rc==GEOMLIB_MESH_SUCCESS && fit.rms<1e-3f &&
          fit.normal.x/fit.normal.z < -0.4f &&
          fit.normal.x/fit.normal.z > -0.6f,"tilted plane fits");
    if(failures!=0){std::printf("%d plane fit test(s) failed\n",failures);return 1;}
    std::printf("All plane fit tests passed\n");
    return 0;
}
