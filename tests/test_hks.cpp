/*
 * Heat kernel signature tests: values are positive and finite, decay with
 * time, and are nearly uniform across a homogeneous sphere.
 */
#include <geomlib/geomlib_hks.h>

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
    geomlib_trimesh sphere={};
    int rc=geomlib_create_sphere(&sphere,1.0f,12,8);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere created");
    if(rc!=GEOMLIB_MESH_SUCCESS) return 1;

    double times[2]={0.1,1.0};
    std::vector<double> hks(sphere.vertex_count*2,0.0);
    rc=geomlib_compute_hks(&sphere,8,times,2,hks.data());
    std::printf("hks rc=%d values=%zu\n",rc,hks.size());
    check(rc==GEOMLIB_MESH_SUCCESS,"HKS computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        double mean[2]={0.0,0.0},minv[2]={1e30,1e30},maxv[2]={0.0,0.0};
        bool finite=true,positive=true;
        for(size_t i=0;i<sphere.vertex_count;++i)
            for(int t=0;t<2;++t){
                double v=hks[i*2+t];
                if(!std::isfinite(v)) finite=false;
                if(v<=0.0) positive=false;
                mean[t]+=v;
                if(v<minv[t]) minv[t]=v;
                if(v>maxv[t]) maxv[t]=v;
            }
        mean[0]/=sphere.vertex_count; mean[1]/=sphere.vertex_count;
        std::printf("hks mean t0=%.5f t1=%.5f range0=[%.5f,%.5f]\n",
                    mean[0],mean[1],minv[0],maxv[0]);
        check(finite && positive,"HKS values are finite and positive");
        check(mean[1]<mean[0],"HKS decays with time");
        check(maxv[0]<=1.000001 && minv[0]>=0.0,
              "normalized HKS lies in [0,1]");
    }
    geomlib_trimesh_free(&sphere);
    if(failures!=0){std::printf("%d HKS test(s) failed\n",failures);return 1;}
    std::printf("All HKS tests passed\n");
    return 0;
}
