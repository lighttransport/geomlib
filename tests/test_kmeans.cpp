/*
 * K-means: two well-separated blobs are recovered with valid labels and
 * centroids near the blob centres, deterministically.
 */
#include <geomlib/geomlib_kmeans.h>

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
    geomlib_pointcloud cloud;
    geomlib_pointcloud_init(&cloud,200);
    for(int i=0;i<200;++i){
        float cx=(i<100)?-3.0f:3.0f;
        float x=cx+0.05f*sinf((float)i*1.7f);
        float y=0.05f*cosf((float)i*2.1f);
        geomlib_pointcloud_add_point(&cloud,geomlib_vec3_create(x,y,0),nullptr,1.0f);
    }
    std::vector<uint32_t> labels(cloud.point_count);
    geomlib_vec3 cent[2];
    int k=geomlib_pointcloud_kmeans(&cloud,2,30,labels.data(),cent);
    size_t c0=0,c1=0;
    for(uint32_t l:labels){ if(l==0)c0++; else if(l==1)c1++; }
    std::printf("kmeans k=%d sizes=%zu,%zu c0=(%.2f,%.2f) c1=(%.2f,%.2f)\n",
                k,c0,c1,cent[0].x,cent[0].y,cent[1].x,cent[1].y);
    check(k==2,"two clusters returned");
    check(c0==100 && c1==100,"clusters split by blob");
    check(std::fabs(cent[0].x+3.0f)<0.2f && std::fabs(cent[1].x-3.0f)<0.2f,
          "centroids match blob centres");

    std::vector<uint32_t> again(cloud.point_count);
    geomlib_pointcloud_kmeans(&cloud,2,30,again.data(),nullptr);
    check(again==labels,"clustering is deterministic");
    geomlib_pointcloud_free(&cloud);
    if(failures!=0){std::printf("%d kmeans test(s) failed\n",failures);return 1;}
    std::printf("All kmeans tests passed\n");
    return 0;
}
