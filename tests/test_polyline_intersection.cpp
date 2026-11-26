/*
 * Polyline intersection: a line through a sphere centre gives two surface
 * crossings; a missing line gives none.
 */
#include <geomlib/geomlib_polyline_intersection.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_trimesh sphere={};
    int rc=geomlib_create_sphere(&sphere,1.0f,16,12);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere created");
    if(rc!=GEOMLIB_MESH_SUCCESS) return 1;

    geomlib_vec3 line[2]={geomlib_vec3_create(-2,0,0),geomlib_vec3_create(2,0,0)};
    geomlib_polyline_hits h=geomlib_polyline_hits_create();
    rc=geomlib_polyline_mesh_intersection(line,2,&sphere,&h);
    std::printf("polyline hit rc=%d count=%zu\n",rc,h.count);
    check(rc==GEOMLIB_MESH_SUCCESS,"intersection computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(h.count==2,"line through centre crosses twice");
        bool on_sphere=true;
        for(size_t i=0;i<h.count;++i)
            if(std::fabs(geomlib_vec3_length(h.points[i])-1.0f)>0.02f)
                on_sphere=false;
        check(on_sphere,"crossings lie on the sphere");
    }
    geomlib_polyline_hits_free(&h);

    geomlib_vec3 miss[2]={geomlib_vec3_create(2,2,2),geomlib_vec3_create(2,2,5)};
    h=geomlib_polyline_hits_create();
    rc=geomlib_polyline_mesh_intersection(miss,2,&sphere,&h);
    std::printf("polyline miss rc=%d count=%zu\n",rc,h.count);
    check(rc==GEOMLIB_MESH_SUCCESS && h.count==0,"missing line gives no hits");
    geomlib_polyline_hits_free(&h);

    /* Three-segment polyline that enters and exits the sphere. */
    geomlib_vec3 poly[3]={geomlib_vec3_create(-2,0,0),geomlib_vec3_create(0,0,0),
                          geomlib_vec3_create(2,0,0)};
    h=geomlib_polyline_hits_create();
    rc=geomlib_polyline_mesh_intersection(poly,3,&sphere,&h);
    check(rc==GEOMLIB_MESH_SUCCESS && h.count==2,"multi-segment polyline works");
    geomlib_polyline_hits_free(&h);

    geomlib_trimesh_free(&sphere);
    if(failures!=0){std::printf("%d polyline intersection test(s) failed\n",failures);return 1;}
    std::printf("All polyline intersection tests passed\n");
    return 0;
}
