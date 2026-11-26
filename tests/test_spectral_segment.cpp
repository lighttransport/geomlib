/*
 * Spectral segmentation: a planar grid splits into two spatial clusters.
 */
#include <geomlib/geomlib_spectral_segment.h>

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
    geomlib_trimesh grid;
    geomlib_trimesh_init(&grid,49,72);
    for(int y=0;y<7;++y)
        for(int x=0;x<7;++x)
            geomlib_trimesh_add_vertex(&grid,geomlib_vec3_create((float)x,(float)y,0),
                                       nullptr,nullptr);
    for(int y=0;y<6;++y)
        for(int x=0;x<6;++x){
            uint32_t i=(uint32_t)(y*7+x);
            geomlib_trimesh_add_triangle(&grid,i,i+1,i+8);
            geomlib_trimesh_add_triangle(&grid,i,i+8,i+7);
        }
    std::vector<uint32_t> labels(grid.vertex_count,0);
    int rc=geomlib_spectral_segment(&grid,2,4,labels.data());
    double sx[2]={0,0};
    size_t cnt[2]={0,0};
    bool valid=true;
    for(size_t i=0;i<grid.vertex_count;++i){
        if(labels[i]>=2) valid=false;
        sx[labels[i]]+=grid.vertices[i].x;
        cnt[labels[i]]++;
    }
    if(cnt[0]) sx[0]/=cnt[0];
    if(cnt[1]) sx[1]/=cnt[1];
    std::printf("spectral rc=%d clusters=%zu,%zu mean x=%.3f,%.3f\n",
                rc,cnt[0],cnt[1],sx[0],sx[1]);
    check(rc==GEOMLIB_MESH_SUCCESS,"spectral segmentation computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(valid,"all labels valid");
        check(cnt[0]>0 && cnt[1]>0,"both clusters non-empty");
        check(std::fabs(sx[0]-sx[1])>1.0,"clusters separate along x");
    }
    geomlib_trimesh_free(&grid);
    if(failures!=0){std::printf("%d spectral test(s) failed\n",failures);return 1;}
    std::printf("All spectral tests passed\n");
    return 0;
}
