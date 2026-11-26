/*
 * 2D Poisson sampling in a polygon: all points inside, separation respected,
 * deterministic.
 */
#include <geomlib/geomlib_poisson2d.h>
#include <geomlib/geomlib_point_in_polygon.h>

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
    geomlib_vec2 box[4]={
        geomlib_vec2_create(0,0),geomlib_vec2_create(2,0),
        geomlib_vec2_create(2,2),geomlib_vec2_create(0,2)
    };
    geomlib_poisson2d_params p;
    p.radius=0.25f; p.candidates=30; p.seed=777; p.max_samples=500;
    std::vector<geomlib_vec2> pts(p.max_samples);
    size_t n=geomlib_poisson_disk_polygon(box,4,&p,pts.data(),pts.size());
    std::printf("poisson2d count=%zu\n",n);
    check(n>20,"samples generated");
    if(n>1){
        float sep=1e9f;
        bool inside=true;
        for(size_t i=0;i<n;++i){
            if(!geomlib_point_in_polygon(box,4,pts[i])) inside=false;
            for(size_t j=0;j<i;++j)
                sep=std::fmin(sep,geomlib_vec2_length(
                    geomlib_vec2_sub(pts[i],pts[j])));
        }
        std::printf("poisson2d min separation=%.4f\n",sep);
        check(inside,"all samples inside polygon");
        check(sep>=p.radius*0.999f,"separation respected");
    }
    std::vector<geomlib_vec2> again(pts.size());
    geomlib_poisson_disk_polygon(box,4,&p,again.data(),again.size());
    bool same=true;
    for(size_t i=0;i<pts.size();++i)
        if(geomlib_vec2_length(geomlib_vec2_sub(pts[i],again[i]))>1e-6f) same=false;
    check(same,"deterministic");
    if(failures!=0){std::printf("%d poisson2d test(s) failed\n",failures);return 1;}
    std::printf("All poisson2d tests passed\n");
    return 0;
}
