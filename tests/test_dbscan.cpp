/*
 * DBSCAN: two dense blobs plus sparse noise yield two clusters and noise
 * labels.
 */
#include <geomlib/geomlib_dbscan.h>

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
    geomlib_pointcloud_init(&cloud,120);
    for(int i=0;i<50;++i)
        geomlib_pointcloud_add_point(&cloud,
            geomlib_vec3_create(-3.0f+0.05f*sinf((float)i),0.05f*cosf((float)i),0),
            nullptr,1.0f);
    for(int i=0;i<50;++i)
        geomlib_pointcloud_add_point(&cloud,
            geomlib_vec3_create(3.0f+0.05f*sinf((float)i),0.05f*cosf((float)i),0),
            nullptr,1.0f);
    for(int i=0;i<20;++i)
        geomlib_pointcloud_add_point(&cloud,
            geomlib_vec3_create((float)(i*7-70),20.0f,0),nullptr,1.0f);

    geomlib_dbscan_params p=geomlib_dbscan_params_default();
    p.eps=0.5f;
    p.min_points=5;
    std::vector<uint32_t> labels(cloud.point_count);
    size_t clusters=0;
    int rc=geomlib_pointcloud_dbscan(&cloud,&p,labels.data(),&clusters);
    size_t noise=0;
    for(uint32_t l:labels) if(l==GEOMLIB_DBSCAN_NOISE) noise++;
    std::printf("dbscan rc=%d clusters=%zu noise=%zu\n",rc,clusters,noise);
    check(rc==GEOMLIB_MESH_SUCCESS,"DBSCAN computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(clusters==2,"two dense blobs found");
        check(noise==20,"sparse points marked as noise");
        bool inner=true;
        for(size_t i=0;i<100;++i)
            if(labels[i]==GEOMLIB_DBSCAN_NOISE) inner=false;
        check(inner,"blob points are clustered");
    }
    geomlib_pointcloud_free(&cloud);
    if(failures!=0){std::printf("%d DBSCAN test(s) failed\n",failures);return 1;}
    std::printf("All DBSCAN tests passed\n");
    return 0;
}
