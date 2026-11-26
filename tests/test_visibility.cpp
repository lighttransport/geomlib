/*
 * Visibility shortest path: a square obstacle in the middle is routed around,
 * while an empty field gives the straight-line distance.
 */
#include <geomlib/geomlib_visibility.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_vec2 s=geomlib_vec2_create(-3,0),g=geomlib_vec2_create(3,0);
    geomlib_visibility_path p=geomlib_visibility_path_create();
    int rc=geomlib_visibility_shortest_path(&s,&g,nullptr,nullptr,0,&p);
    std::printf("empty rc=%d len=%.4f nodes=%zu\n",rc,p.length,p.count);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(p.length-6.0)<1e-4,
          "empty field gives straight distance");
    geomlib_visibility_path_free(&p);

    geomlib_vec2 box[4]={
        geomlib_vec2_create(-1,-1),geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1),geomlib_vec2_create(-1,1)
    };
    const geomlib_vec2 *obs[1]={box};
    size_t sizes[1]={4};
    rc=geomlib_visibility_shortest_path(&s,&g,obs,sizes,1,&p);
    std::printf("obstacle rc=%d len=%.4f nodes=%zu\n",rc,p.length,p.count);
    check(rc==GEOMLIB_MESH_SUCCESS,"path found around obstacle");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(p.length>6.0f && p.length<9.0f,"path detours around the square");
        check(p.count>=3,"path bends at obstacle corners");
        /* All segments must avoid the obstacle interior midpoint check. */
        geomlib_visibility_path_free(&p);
    }
    if(failures!=0){std::printf("%d visibility test(s) failed\n",failures);return 1;}
    std::printf("All visibility tests passed\n");
    return 0;
}
