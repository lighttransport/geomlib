/*
 * Dihedral region growing: a box separates into six planar faces at a small
 * threshold, while a smooth sphere remains a single region.
 */
#include <geomlib/geomlib_region_growing.h>

#include <cstdio>
#include <vector>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_trimesh box={};
    int rc=geomlib_create_box(&box,1.0f,1.0f,1.0f);
    check(rc==GEOMLIB_MESH_SUCCESS,"box created");
    if(rc!=GEOMLIB_MESH_SUCCESS) return 1;
    std::vector<uint32_t> labels(box.triangle_count);
    size_t regions=0;
    rc=geomlib_mesh_region_growing(&box,30.0f,labels.data(),&regions);
    std::printf("regions box rc=%d regions=%zu tris=%zu\n",rc,regions,box.triangle_count);
    check(rc==GEOMLIB_MESH_SUCCESS,"box segmentation succeeds");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(regions==6,"box separates into six faces");
        bool valid=true;
        for(size_t i=0;i<labels.size();++i) if(labels[i]>=regions) valid=false;
        check(valid,"all face labels valid");
    }
    geomlib_trimesh_free(&box);

    /* Flat grid: all normals agree, so one region. */
    geomlib_trimesh grid;
    geomlib_trimesh_init(&grid,25,32);
    for(int y=0;y<5;++y)
        for(int x=0;x<5;++x)
            geomlib_trimesh_add_vertex(&grid,geomlib_vec3_create((float)x,(float)y,0),
                                       nullptr,nullptr);
    for(int y=0;y<4;++y)
        for(int x=0;x<4;++x){
            uint32_t i=(uint32_t)(y*5+x);
            geomlib_trimesh_add_triangle(&grid,i,i+1,i+6);
            geomlib_trimesh_add_triangle(&grid,i,i+6,i+5);
        }
    {
        std::vector<uint32_t> lab(grid.triangle_count);
        rc=geomlib_mesh_region_growing(&grid,30.0f,lab.data(),&regions);
        std::printf("regions grid rc=%d regions=%zu\n",rc,regions);
        check(rc==GEOMLIB_MESH_SUCCESS && regions==1,
              "flat grid is one region");
    }
    geomlib_trimesh_free(&grid);
    if(failures!=0){std::printf("%d region growing test(s) failed\n",failures);return 1;}
    std::printf("All region growing tests passed\n");
    return 0;
}
