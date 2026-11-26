/*
 * 3D Poisson disk: all points in box, pairwise separation >= radius and
 * deterministic for a fixed seed.
 */
#include <geomlib/geomlib_poisson3d.h>

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
    geomlib_poisson3d_params p;
    p.box_min=geomlib_vec3_create(0,0,0);
    p.box_max=geomlib_vec3_create(2,2,2);
    p.radius=0.25f;
    p.candidates=30;
    p.seed=12345;
    p.max_samples=2000;
    std::vector<geomlib_vec3> pts(p.max_samples);
    size_t n=geomlib_poisson_disk_volume(&p,pts.data(),pts.size());
    std::printf("poisson3d count=%zu\n",n);
    check(n>20,"samples generated");
    if(n>1){
        float minsep=1e9f;
        bool in_box=true;
        for(size_t i=0;i<n;++i){
            geomlib_vec3 q=pts[i];
            if(q.x<0||q.x>2||q.y<0||q.y>2||q.z<0||q.z>2) in_box=false;
            for(size_t j=0;j<i;++j)
                minsep=std::fmin(minsep,geomlib_vec3_length(geomlib_vec3_sub(pts[i],pts[j])));
        }
        std::printf("poisson3d min separation=%.4f\n",minsep);
        check(in_box,"all points inside the box");
        check(minsep>=p.radius*0.999f,"minimum separation respected");
    }
    std::vector<geomlib_vec3> again(pts.size());
    geomlib_poisson_disk_volume(&p,again.data(),again.size());
    bool same=true;
    for(size_t i=0;i<pts.size();++i)
        if(geomlib_vec3_length(geomlib_vec3_sub(pts[i],again[i]))>1e-6f) same=false;
    check(same,"deterministic for fixed seed");
    if(failures!=0){std::printf("%d poisson3d test(s) failed\n",failures);return 1;}
    std::printf("All poisson3d tests passed\n");
    return 0;
}
