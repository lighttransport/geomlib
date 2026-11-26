/*
 * Geodesic diameter: a unit sphere's graph diameter is near pi; the exact
 * and double-sweep results agree and endpoints are valid.
 */
#include <geomlib/geomlib_geodesic_diameter.h>
#include <geomlib/geomlib_geodesic_dijkstra.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_trimesh m={};
    int rc=geomlib_create_sphere(&m,1.0f,10,7);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere created");
    if(rc!=GEOMLIB_MESH_SUCCESS) return 1;

    geomlib_geodesic_diameter exact,sweep;
    rc=geomlib_geodesic_diameter_exact(&m,&exact);
    check(rc==GEOMLIB_MESH_SUCCESS,"exact diameter computes");
    rc=geomlib_geodesic_diameter_sweep(&m,&sweep);
    check(rc==GEOMLIB_MESH_SUCCESS,"sweep diameter computes");
    std::printf("diameter exact=%.4f sweep=%.4f endpoints=%u,%u\n",
                exact.distance,sweep.distance,exact.source,exact.target);
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(exact.source<m.vertex_count && exact.target<m.vertex_count,
              "endpoints are valid");
        check(exact.distance>2.5f && exact.distance<3.6f,
              "sphere diameter is near pi");
        check(std::fabs(sweep.distance-exact.distance)<0.2f,
              "double sweep agrees with exact");
        geomlib_geodesic_distances d=geomlib_geodesic_create(m.vertex_count);
        geomlib_geodesic_dijkstra_from_vertex(&m,exact.source,&d);
        check(std::fabs(d.distances[exact.target]-exact.distance)<1e-3f,
              "reported distance matches Dijkstra");
        geomlib_geodesic_free(&d);
    }
    geomlib_trimesh_free(&m);
    if(failures!=0){std::printf("%d diameter test(s) failed\n",failures);return 1;}
    std::printf("All diameter tests passed\n");
    return 0;
}
