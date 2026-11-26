/*
 * Mesh-mesh distance: two separated spheres give the analytic gap, and
 * overlapping meshes give zero.
 */
#include <geomlib/geomlib_mesh_distance.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_trimesh a={},b={};
    int rc=geomlib_create_sphere(&a,1.0f,16,12);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere A");
    rc=geomlib_create_sphere(&b,1.0f,16,12);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere B");
    for(size_t i=0;i<b.vertex_count;++i) b.vertices[i].x+=3.0f;

    geomlib_mesh_distance d;
    rc=geomlib_mesh_mesh_distance(&a,&b,&d);
    std::printf("distance separated rc=%d d=%.4f a=(%.3f,%.3f,%.3f) b=(%.3f,%.3f,%.3f)\n",rc,d.distance,d.point_a.x,d.point_a.y,d.point_a.z,d.point_b.x,d.point_b.y,d.point_b.z);
    check(rc==GEOMLIB_MESH_SUCCESS,"distance computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(d.distance-1.0f)<0.05f,"analytic gap is 1");
        check(std::fabs(d.point_b.x-d.point_a.x-1.0f)<0.05f &&
              std::fabs(d.point_a.y)<0.05f && std::fabs(d.point_a.z)<0.05f,
              "closest points face each other along x");
    }

    for(size_t i=0;i<b.vertex_count;++i) b.vertices[i].x-=1.5f;
    rc=geomlib_mesh_mesh_distance(&a,&b,&d);
    std::printf("distance overlapping rc=%d d=%.4f\n",rc,d.distance);
    check(rc==GEOMLIB_MESH_SUCCESS && d.distance<0.05f,
          "overlapping meshes have zero distance");

    geomlib_trimesh_free(&a);
    geomlib_trimesh_free(&b);
    if(failures!=0){std::printf("%d mesh distance test(s) failed\n",failures);return 1;}
    std::printf("All mesh distance tests passed\n");
    return 0;
}
