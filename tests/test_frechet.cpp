/*
 * Discrete Fréchet: identical curves give 0, translated curves give the
 * translation, and the distance is symmetric.
 */
#include <geomlib/geomlib_frechet.h>

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
    std::vector<geomlib_vec3> a={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(2,0,0)
    };
    double d=geomlib_frechet_distance(a.data(),a.size(),a.data(),a.size());
    std::printf("frechet identical=%.6f\n",d);
    check(d<1e-6,"identical curves give zero");

    std::vector<geomlib_vec3> b;
    for(auto &p:a) b.push_back(geomlib_vec3_create(p.x,p.y+3.0f,p.z));
    d=geomlib_frechet_distance(a.data(),a.size(),b.data(),b.size());
    std::printf("frechet translated=%.6f\n",d);
    check(std::fabs(d-3.0)<1e-5,"translation distance is 3");
    double back=geomlib_frechet_distance(b.data(),b.size(),a.data(),a.size());
    check(std::fabs(d-back)<1e-6,"symmetric");
    if(failures!=0){std::printf("%d Fréchet test(s) failed\n",failures);return 1;}
    std::printf("All Fréchet tests passed\n");
    return 0;
}
