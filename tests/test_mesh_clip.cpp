/*
 * Plane clipping: half a box gives half the volume and a capped closed mesh.
 */
#include <geomlib/geomlib_mesh_clip.h>
#include <geomlib/geomlib_mass_properties.h>
#include <geomlib/geomlib_topology.h>

#include <cmath>
#include <cstdio>

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

    geomlib_trimesh out;
    rc=geomlib_clip_mesh_plane(&box,geomlib_vec3_create(0,0,0),
                               geomlib_vec3_create(1,0,0),0,&out);
    std::printf("clip rc=%d verts=%zu tris=%zu\n",rc,out.vertex_count,out.triangle_count);
    check(rc==GEOMLIB_MESH_SUCCESS,"clipping succeeds");
    if(rc==GEOMLIB_MESH_SUCCESS){
        geomlib_mass_properties mp;
        rc=geomlib_compute_mass_properties(&out,&mp);
        std::printf("clip volume=%.4f expect=0.5\n",mp.volume);
        check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(mp.volume-0.5)<0.02,
              "half box volume is 0.5");
        geomlib_trimesh_free(&out);
    }
    geomlib_trimesh_free(&box);
    if(failures!=0){std::printf("%d mesh clip test(s) failed\n",failures);return 1;}
    std::printf("All mesh clip tests passed\n");
    return 0;
}
