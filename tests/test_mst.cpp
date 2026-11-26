/*
 * MST: n-1 edges, connectivity, and known total weight for a unit square.
 */
#include <geomlib/geomlib_mst.h>

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
    double sq[8]={0,0, 1,0, 1,1, 0,1};
    geomlib_mst_result r=geomlib_mst_result_create();
    int rc=geomlib_compute_mst(sq,4,2,&r);
    std::printf("mst rc=%d edges=%zu weight=%.4f expect=3\n",
                rc,r.edge_count,r.total_weight);
    check(rc==GEOMLIB_MESH_SUCCESS && r.edge_count==3 &&
          std::fabs(r.total_weight-3.0)<1e-5,"square MST weight is 3");

    /* Connectedness via union-find. */
    int parent[4]={0,1,2,3};
    auto find=[&](int x){while(parent[x]!=x){parent[x]=parent[parent[x]];x=parent[x];}return x;};
    for(size_t i=0;i<r.edge_count;++i){
        int a=find(r.edges[i*2]),b=find(r.edges[i*2+1]);
        if(a!=b) parent[a]=b;
    }
    bool connected=true;
    for(int i=1;i<4;++i) if(find(i)!=find(0)) connected=false;
    check(connected,"MST connects all vertices");
    geomlib_mst_result_free(&r);

    double line[6]={0,0, 1,0, 3,0};
    rc=geomlib_compute_mst(line,3,2,&r);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(r.total_weight-3.0)<1e-5,
          "collinear MST weight is 3");
    geomlib_mst_result_free(&r);
    if(failures!=0){std::printf("%d MST test(s) failed\n",failures);return 1;}
    std::printf("All MST tests passed\n");
    return 0;
}
