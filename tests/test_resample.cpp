/*
 * Arc-length resampling: endpoints preserved, spacing uniform, and total
 * polyline length recovered.
 */
#include <geomlib/geomlib_resample.h>

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
    std::vector<geomlib_vec3> poly={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(3,0,0),
        geomlib_vec3_create(3.5f,0,0),geomlib_vec3_create(9,0,0)
    };
    std::vector<geomlib_vec3> out(21);
    int rc=geomlib_resample_polyline(poly.data(),poly.size(),out.size(),out.data());
    check(rc==GEOMLIB_MESH_SUCCESS,"resampling succeeds");
    if(rc==GEOMLIB_MESH_SUCCESS){
        double mn=1e30,mx=0,total=0.0;
        for(size_t i=1;i<out.size();++i){
            float d=geomlib_vec3_length(geomlib_vec3_sub(out[i],out[i-1]));
            total+=d;
            if(d<mn)mn=d; if(d>mx)mx=d;
        }
        std::printf("resample count=%zu length=%.4f expect=9 spacing=[%.4f,%.4f]\n",
                    out.size(),total,mn,mx);
        check(std::fabs(total-9.0)<1e-3,"total length preserved");
        check(mx-mn<1e-3,"spacing is uniform");
        check(geomlib_vec3_length(geomlib_vec3_sub(out[0],poly[0]))<1e-5 &&
              geomlib_vec3_length(geomlib_vec3_sub(out.back(),poly.back()))<1e-5,
              "endpoints preserved");
    }
    if(failures!=0){std::printf("%d resample test(s) failed\n",failures);return 1;}
    std::printf("All resample tests passed\n");
    return 0;
}
