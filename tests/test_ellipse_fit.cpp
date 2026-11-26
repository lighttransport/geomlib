/*
 * Direct ellipse fit: recovers centre, axes and angle of a rotated ellipse.
 */
#include <geomlib/geomlib_ellipse_fit.h>

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
    const float cx=1.0f,cy=-2.0f,A=3.0f,B=1.5f,ang=0.4f;
    std::vector<geomlib_vec2> pts;
    for(int i=0;i<200;++i){
        float t=i*6.283185f/200.0f;
        float u=A*std::cos(t),v=B*std::sin(t);
        pts.push_back(geomlib_vec2_create(
            cx+u*std::cos(ang)-v*std::sin(ang),
            cy+u*std::sin(ang)+v*std::cos(ang)));
    }
    geomlib_ellipse_fit fit;
    int rc=geomlib_fit_ellipse(pts.data(),pts.size(),&fit);
    std::printf("fit ellipse rc=%d c=(%.3f,%.3f) axes=(%.3f,%.3f) angle=%.3f rms=%.5f\n",
                rc,fit.center.x,fit.center.y,fit.axis_a,fit.axis_b,fit.angle,fit.rms);
    check(rc==GEOMLIB_MESH_SUCCESS,"ellipse fit computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(fit.center.x-cx)<0.05f && std::fabs(fit.center.y-cy)<0.05f,
              "centre recovered");
        check(std::fabs(fit.axis_a-A)<0.05f && std::fabs(fit.axis_b-B)<0.05f,
              "semi-axes recovered");
        check(fit.rms<0.02f,"low residual");
    }
    if(failures!=0){std::printf("%d ellipse fit test(s) failed\n",failures);return 1;}
    std::printf("All ellipse fit tests passed\n");
    return 0;
}
