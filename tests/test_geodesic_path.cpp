/*
 * Geodesic path: endpoints correct, consecutive vertices are mesh edges, and
 * path length matches the Dijkstra distance.
 */
#include <geomlib/geomlib_geodesic_path.h>
#include <geomlib/geomlib_geodesic_dijkstra.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static bool is_edge(const geomlib_trimesh &m, uint32_t a, uint32_t b)
{
    for(size_t i=0;i<m.triangle_count;++i){
        uint32_t v[3]={m.indices[i*3],m.indices[i*3+1],m.indices[i*3+2]};
        for(int k=0;k<3;++k)
            if((v[k]==a&&v[(k+1)%3]==b)||(v[k]==b&&v[(k+1)%3]==a)) return true;
    }
    return false;
}

int main(void)
{
    geomlib_trimesh m={};
    int rc=geomlib_create_sphere(&m,1.0f,12,8);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere created");
    if(rc!=GEOMLIB_MESH_SUCCESS) return 1;

    uint32_t s=0,t=(uint32_t)(m.vertex_count/2);
    geomlib_geodesic_path p=geomlib_geodesic_path_create();
    rc=geomlib_geodesic_path_between(&m,s,t,&p);
    std::printf("path rc=%d vertices=%zu length=%.4f\n",rc,p.count,p.length);
    check(rc==GEOMLIB_MESH_SUCCESS,"path computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(p.count>=2 && p.vertices[0]==s && p.vertices[p.count-1]==t,
              "endpoints are source and target");
        bool edges=true;
        for(size_t i=1;i<p.count;++i)
            if(!is_edge(m,p.vertices[i-1],p.vertices[i])) edges=false;
        check(edges,"consecutive vertices are mesh edges");
        geomlib_geodesic_distances d=geomlib_geodesic_create(m.vertex_count);
        geomlib_geodesic_dijkstra_from_vertex(&m,s,&d);
        std::printf("path dijkstra=%.4f\n",d.distances[t]);
        check(std::fabs(d.distances[t]-p.length)<1e-3f,
              "path length matches Dijkstra distance");
        geomlib_geodesic_free(&d);
        geomlib_geodesic_path_free(&p);
    }

    rc=geomlib_geodesic_path_between(&m,0,0,&p);
    check(rc==GEOMLIB_MESH_SUCCESS && p.count==1 && p.length==0.0f,
          "source equals target gives trivial path");
    geomlib_geodesic_path_free(&p);
    geomlib_trimesh_free(&m);
    if(failures!=0){std::printf("%d geodesic path test(s) failed\n",failures);return 1;}
    std::printf("All geodesic path tests passed\n");
    return 0;
}
