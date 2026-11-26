/*
 * Edge-collapse coarsening: removes short edges on a refined sphere while
 * keeping the mesh manifold and its boundary unchanged.
 */
#include <geomlib/geomlib_remesh_split.h>
#include <geomlib/geomlib_remesh_collapse.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static bool manifold(const geomlib_trimesh &m)
{
    std::vector<uint64_t> e;
    for(size_t i=0;i<m.triangle_count;++i)
        for(int k=0;k<3;++k){
            uint32_t a=m.indices[i*3+k],b=m.indices[i*3+(k+1)%3];
            if(a>b) std::swap(a,b);
            e.push_back(((uint64_t)a<<32)|b);
        }
    std::sort(e.begin(),e.end());
    for(size_t i=0;i<e.size();){
        size_t j=i+1; while(j<e.size()&&e[j]==e[i]) ++j;
        if(j-i>2) return false;
        i=j;
    }
    return true;
}

int main(void)
{
    geomlib_trimesh m={};
    int rc=geomlib_create_sphere(&m,1.0f,12,8);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere created");
    if(rc!=GEOMLIB_MESH_SUCCESS) return 1;
    geomlib_refine_edges(&m,0.30f,400);
    size_t t0=m.triangle_count,v0=m.vertex_count;
    float vol0=geomlib_trimesh_volume(&m);

    int collapses=geomlib_coarsen_edges(&m,0.2f,400,1);
    std::printf("collapse count=%d verts %zu->%zu tris %zu->%zu vol %.4f->%.4f\n",
                collapses,v0,m.vertex_count,t0,m.triangle_count,vol0,
                geomlib_trimesh_volume(&m));
    check(collapses>0,"collapses performed");
    check(m.triangle_count<t0 && m.vertex_count<v0,"mesh was coarsened");
    check(manifold(m),"mesh remains manifold");
    check(std::fabs(geomlib_trimesh_volume(&m)-vol0)<0.3f,
          "volume approximately preserved");
    geomlib_trimesh_free(&m);
    if(failures!=0){std::printf("%d collapse test(s) failed\n",failures);return 1;}
    std::printf("All collapse tests passed\n");
    return 0;
}
