/*
 * Geodesic ball: radius zero includes only the source, larger radii grow
 * monotonically and never exceed the vertex count.
 */
#include <geomlib/geomlib_geodesic_ball.h>
#include <geomlib/geomlib_geodesic_dijkstra.h>

#include <cstdio>
#include <vector>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_trimesh m;
    geomlib_trimesh_init(&m,25,32);
    for(int y=0;y<5;++y)
        for(int x=0;x<5;++x)
            geomlib_trimesh_add_vertex(&m,geomlib_vec3_create((float)x,(float)y,0),
                                       nullptr,nullptr);
    for(int y=0;y<4;++y)
        for(int x=0;x<4;++x){
            uint32_t i=(uint32_t)(y*5+x);
            geomlib_trimesh_add_triangle(&m,i,i+1,i+6);
            geomlib_trimesh_add_triangle(&m,i,i+6,i+5);
        }
    geomlib_geodesic_ball b=geomlib_geodesic_ball_create();
    int rc=geomlib_extract_geodesic_ball(&m,0,0.0f,&b);
    std::printf("ball zero rc=%d count=%zu\n",rc,b.count);
    check(rc==GEOMLIB_MESH_SUCCESS && b.count==1 &&
          b.vertices[0]==0 && b.mask[0]==1,"zero radius is the source only");
    geomlib_geodesic_ball_free(&b);

    size_t counts[3];
    float radii[3]={0.5f,1.0f,2.0f};
    for(int i=0;i<3;++i){
        b=geomlib_geodesic_ball_create();
        rc=geomlib_extract_geodesic_ball(&m,0,radii[i],&b);
        counts[i]=b.count;
        check(rc==GEOMLIB_MESH_SUCCESS && b.count>0 && b.count<=m.vertex_count,
              "radius query is valid");
        bool consistent=true;
        for(size_t v=0;v<m.vertex_count;++v)
            if(b.mask[v] && b.count==0) consistent=false;
        check(consistent,"mask allocated");
        geomlib_geodesic_ball_free(&b);
    }
    std::printf("ball counts r=0.5,1.0,2.0: %zu,%zu,%zu of %zu\n",
                counts[0],counts[1],counts[2],m.vertex_count);
    check(counts[0]<=counts[1] && counts[1]<=counts[2],"ball grows with radius");

    b=geomlib_geodesic_ball_create();
    rc=geomlib_extract_geodesic_ball(&m,0,10.0f,&b);
    check(rc==GEOMLIB_MESH_SUCCESS && b.count==m.vertex_count,
          "large radius includes all vertices");
    geomlib_geodesic_ball_free(&b);

    geomlib_trimesh_free(&m);
    if(failures!=0){std::printf("%d geodesic ball test(s) failed\n",failures);return 1;}
    std::printf("All geodesic ball tests passed\n");
    return 0;
}
