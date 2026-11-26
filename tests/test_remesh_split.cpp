/*
 * Edge-split refinement: bounding the longest edge while preserving volume
 * and manifoldness of a sphere.
 */
#include <geomlib/geomlib_remesh_split.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static float max_edge(const geomlib_trimesh &m)
{
    float mx=0.0f;
    for(size_t i=0;i<m.triangle_count;++i)
        for(int k=0;k<3;++k){
            uint32_t a=m.indices[i*3+k],b=m.indices[i*3+(k+1)%3];
            mx=std::fmax(mx,geomlib_vec3_length(geomlib_vec3_sub(m.vertices[a],m.vertices[b])));
        }
    return mx;
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
    int rc=geomlib_create_sphere(&m,1.0f,8,6);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere created");
    if(rc!=GEOMLIB_MESH_SUCCESS) return 1;
    float before=max_edge(m);
    size_t t0=m.triangle_count,v0=m.vertex_count;
    float vol0=geomlib_trimesh_volume(&m);

    int splits=geomlib_refine_edges(&m,0.35f,400);
    float after=max_edge(m);
    float vol1=geomlib_trimesh_volume(&m);
    std::printf("refine splits=%d verts %zu->%zu tris %zu->%zu max edge %.4f->%.4f vol %.4f->%.4f\n",
                splits,v0,m.vertex_count,t0,m.triangle_count,before,after,vol0,vol1);
    check(splits>0,"splits performed");
    check(after<=0.36f,"longest edge is bounded");
    check(m.triangle_count>t0 && m.vertex_count>v0,"mesh was refined");
    check(std::fabs(vol1-vol0)<0.25f,"volume approximately preserved");
    check(manifold(m),"mesh remains manifold");
    geomlib_trimesh_free(&m);
    if(failures!=0){std::printf("%d refine test(s) failed\n",failures);return 1;}
    std::printf("All refine tests passed\n");
    return 0;
}
