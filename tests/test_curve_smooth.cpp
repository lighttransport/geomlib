/*
 * Curve smoothing: endpoints fixed, a noisy zigzag becomes smoother and a
 * straight line stays straight.
 */
#include <geomlib/geomlib_curve_smooth.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static double roughness(const std::vector<geomlib_vec3> &p)
{
    double s=0.0;
    for(size_t i=1;i+1<p.size();++i) {
        geomlib_vec3 d=geomlib_vec3_sub(p[i],geomlib_vec3_mul(
            geomlib_vec3_add(p[i-1],p[i+1]),0.5f));
        s+=geomlib_vec3_length_squared(d);
    }
    return s;
}

int main(void)
{
    std::vector<geomlib_vec3> pts;
    for(int i=0;i<21;++i)
        pts.push_back(geomlib_vec3_create((float)i,
            (i%2)?0.4f:0.0f,0));
    double before=roughness(pts);
    int rc=geomlib_smooth_polyline(pts.data(),pts.size(),100,0.5f);
    double after=roughness(pts);
    std::printf("curve smooth rc=%d roughness %.4f -> %.6f\n",rc,before,after);
    check(rc==GEOMLIB_MESH_SUCCESS,"smoothing succeeds");
    check(after<before*0.01f,"roughness drops");
    check(std::fabs(pts[0].x)<1e-6 && std::fabs(pts.back().x-20.0f)<1e-6,
          "endpoints preserved");

    std::vector<geomlib_vec3> line;
    for(int i=0;i<11;++i) line.push_back(geomlib_vec3_create((float)i,0,0));
    rc=geomlib_smooth_polyline(line.data(),line.size(),10,0.5f);
    bool straight=true;
    for(auto &p:line) if(std::fabs(p.y)>1e-6f) straight=false;
    check(rc==GEOMLIB_MESH_SUCCESS && straight,"straight line stays straight");
    if(failures!=0){std::printf("%d curve smooth test(s) failed\n",failures);return 1;}
    std::printf("All curve smooth tests passed\n");
    return 0;
}
