/*
 * Edge-flip quality improvement: raises the minimum angle on a distorted
 * planar mesh without changing topology, triangle count or boundary.
 */
#include <geomlib/geomlib_quality_flip.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static double area(const geomlib_trimesh &m)
{
    double s=0.0;
    for(size_t i=0;i<m.triangle_count;++i){
        geomlib_vec3 a=m.vertices[m.indices[i*3]];
        geomlib_vec3 b=m.vertices[m.indices[i*3+1]];
        geomlib_vec3 c=m.vertices[m.indices[i*3+2]];
        s+=0.5*geomlib_vec3_length(geomlib_vec3_cross(
            geomlib_vec3_sub(b,a),geomlib_vec3_sub(c,a)));
    }
    return s;
}

static std::vector<uint64_t> boundary_edges(const geomlib_trimesh &m)
{
    std::vector<uint64_t> e;
    for(size_t i=0;i<m.triangle_count;++i){
        uint32_t v[3]={m.indices[i*3],m.indices[i*3+1],m.indices[i*3+2]};
        for(int k=0;k<3;++k){
            uint32_t a=v[k],b=v[(k+1)%3];
            if(a>b) std::swap(a,b);
            e.push_back(((uint64_t)a<<32)|b);
        }
    }
    std::sort(e.begin(),e.end());
    std::vector<uint64_t> out;
    for(size_t i=0;i<e.size();){
        size_t j=i+1;
        while(j<e.size()&&e[j]==e[i]) ++j;
        if(j-i==1) out.push_back(e[i]);
        i=j;
    }
    return out;
}

int main(void)
{
    /* Irregular planar triangulation of a 5x5 grid. */
    geomlib_trimesh m;
    geomlib_trimesh_init(&m,25,32);
    for(int y=0;y<5;++y)
        for(int x=0;x<5;++x){
            float xx=(float)x+((x>0&&x<4&&y>0&&y<4)?0.45f*sinf((float)(x*3+y)):0.0f);
            geomlib_trimesh_add_vertex(&m,geomlib_vec3_create(xx,(float)y+((x>0&&x<4&&y>0&&y<4)?0.3f*cosf((float)(x*5+y*2)):0.0f),0),nullptr,nullptr);
        }
    for(int y=0;y<4;++y)
        for(int x=0;x<4;++x){
            uint32_t i=(uint32_t)(y*5+x);
            if((x+y)%2){geomlib_trimesh_add_triangle(&m,i,i+1,i+6);
                        geomlib_trimesh_add_triangle(&m,i,i+6,i+5);}
            else{geomlib_trimesh_add_triangle(&m,i,i+1,i+5);
                 geomlib_trimesh_add_triangle(&m,i+1,i+6,i+5);}
        }
    float before=geomlib_mesh_min_angle(&m);
    double area0=area(m);
    auto bd0=boundary_edges(m);
    int flips=geomlib_improve_mesh_quality(&m,30);
    float after=geomlib_mesh_min_angle(&m);
    double area1=area(m);
    auto bd1=boundary_edges(m);
    std::printf("quality flips=%d min angle before=%.4f after=%.4f area %.4f->%.4f\n",
                flips,before,after,area0,area1);
    check(flips>0,"edge flips performed");
    check(after>before,"minimum angle improved");
    check(std::fabs(area1-area0)<1e-3,"area preserved");
    std::printf("boundary edges before=%zu after=%zu\n",bd0.size(),bd1.size());
    if (bd0!=bd1) {
        for (auto x:bd0) if (std::find(bd1.begin(),bd1.end(),x)==bd1.end())
            std::printf("  lost %u-%u\n",(uint32_t)(x>>32),(uint32_t)x);
        for (auto x:bd1) if (std::find(bd0.begin(),bd0.end(),x)==bd0.end())
            std::printf("  gained %u-%u\n",(uint32_t)(x>>32),(uint32_t)x);
    }
    check(bd0==bd1,"boundary edges unchanged");
    check(m.triangle_count==32,"triangle count unchanged");
    geomlib_trimesh_free(&m);
    if(failures!=0){std::printf("%d quality flip test(s) failed\n",failures);return 1;}
    std::printf("All quality flip tests passed\n");
    return 0;
}
