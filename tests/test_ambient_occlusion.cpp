/*
 * Ambient occlusion: a flat grid is fully open, a concave box interior has
 * lower visibility, and all values stay in [0,1].
 */
#include <geomlib/geomlib_ambient_occlusion.h>

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
    std::vector<float> ao(grid.vertex_count);
    int rc=geomlib_mesh_ambient_occlusion(&grid,32,10.0f,ao.data());
    float mn=1.0f,mx=0.0f;
    for(float v:ao){ if(v<mn)mn=v; if(v>mx)mx=v; }
    std::printf("ao plane rc=%d min=%.4f max=%.4f\n",rc,mn,mx);
    check(rc==GEOMLIB_MESH_SUCCESS,"plane AO computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(mn>0.95f,"flat plane is fully open");
        check(mx<=1.0f,"AO values bounded by 1");
    }
    geomlib_trimesh_free(&grid);

    /* Inside corner of three perpendicular quads: the origin vertex sees
     * the other two walls, so its AO is below 1. */
    geomlib_trimesh room;
    geomlib_trimesh_init(&room,7,6);
    geomlib_vec3 rv[7]={
        geomlib_vec3_create(0,0,0), geomlib_vec3_create(0,1,0),
        geomlib_vec3_create(0,1,1), geomlib_vec3_create(0,0,1),
        geomlib_vec3_create(1,0,0), geomlib_vec3_create(1,0,1),
        geomlib_vec3_create(1,1,0)
    };
    for(int i=0;i<7;++i) geomlib_trimesh_add_vertex(&room,rv[i],nullptr,nullptr);
    geomlib_trimesh_add_triangle(&room,0,1,2);
    geomlib_trimesh_add_triangle(&room,0,2,3);
    geomlib_trimesh_add_triangle(&room,0,4,5);
    geomlib_trimesh_add_triangle(&room,0,5,3);
    geomlib_trimesh_add_triangle(&room,0,6,4);
    geomlib_trimesh_add_triangle(&room,0,1,6);
    std::vector<float> rao(room.vertex_count);
    rc=geomlib_mesh_ambient_occlusion(&room,64,10.0f,rao.data());
    std::printf("ao room rc=%d origin=%.4f\n",rc,rao[0]);
    check(rc==GEOMLIB_MESH_SUCCESS && rao[0]<0.99f && rao[0]>=0.0f,
          "inside corner is occluded");
    geomlib_trimesh_free(&room);
    if(failures!=0){std::printf("%d AO test(s) failed\n",failures);return 1;}
    std::printf("All AO tests passed\n");
    return 0;
}
