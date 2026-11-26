/*
 * Surface Voronoi tests: every vertex is assigned, sites own themselves,
 * regions are connected and partition the mesh.
 */
#include <geomlib/geomlib_surface_voronoi.h>

#include <cstdio>
#include <vector>
#include <queue>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_trimesh m={};
    int rc=geomlib_create_sphere(&m,1.0f,12,8);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere created");
    if(rc!=GEOMLIB_MESH_SUCCESS) return 1;
    uint32_t sites[3]={0, (uint32_t)(m.vertex_count/3),
                       (uint32_t)(2*m.vertex_count/3)};
    geomlib_surface_voronoi v=geomlib_surface_voronoi_create();
    rc=geomlib_surface_voronoi_compute(&m,sites,3,&v);
    std::printf("voronoi rc=%d vertices=%zu sizes=%zu,%zu,%zu\n",rc,
                v.vertex_count,v.sizes[0],v.sizes[1],v.sizes[2]);
    check(rc==GEOMLIB_MESH_SUCCESS,"surface Voronoi computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        bool assigned=true;
        size_t total=0;
        for(size_t i=0;i<v.vertex_count;++i)
            if(v.labels[i]>=3) assigned=false;
        for(size_t k=0;k<3;++k) total+=v.sizes[k];
        check(assigned,"all labels valid");
        check(total==v.vertex_count,"regions partition the vertices");
        for(size_t k=0;k<3;++k)
            check(v.labels[sites[k]]==(uint32_t)k,"each site owns itself");
        /* Region connectivity via mesh edges. */
        std::vector<std::vector<uint32_t>> adj(m.vertex_count);
        for(size_t i=0;i<m.triangle_count;++i){
            uint32_t a=m.indices[i*3],b=m.indices[i*3+1],c=m.indices[i*3+2];
            adj[a].push_back(b); adj[a].push_back(c);
            adj[b].push_back(a); adj[b].push_back(c);
            adj[c].push_back(a); adj[c].push_back(b);
        }
        bool connected=true;
        for(size_t k=0;k<3;++k){
            std::vector<uint8_t> seen(m.vertex_count,0);
            std::queue<uint32_t> q;
            q.push(sites[k]); seen[sites[k]]=1;
            size_t count=0;
            while(!q.empty()){
                uint32_t u=q.front(); q.pop(); count++;
                for(uint32_t w:adj[u])
                    if(!seen[w]&&v.labels[w]==(uint32_t)k){
                        seen[w]=1; q.push(w);
                    }
            }
            if(count!=v.sizes[k]) connected=false;
        }
        check(connected,"each Voronoi region is connected");
        geomlib_surface_voronoi_free(&v);
    }
    geomlib_trimesh_free(&m);
    if(failures!=0){std::printf("%d surface Voronoi test(s) failed\n",failures);return 1;}
    std::printf("All surface Voronoi tests passed\n");
    return 0;
}
