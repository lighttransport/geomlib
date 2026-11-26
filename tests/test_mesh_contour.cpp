/*
 * Contour tests: a planar grid with a linear field gives straight parallel
 * isolines whose endpoints carry the iso value.
 */
#include <geomlib/geomlib_mesh_contour.h>

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
    std::vector<float> vals(m.vertex_count);
    for(size_t i=0;i<m.vertex_count;++i) vals[i]=m.vertices[i].x;

    geomlib_mesh_contour c=geomlib_mesh_contour_create();
    int rc=geomlib_mesh_contour_lines(&m,vals.data(),2.0f,&c);
    std::printf("contour rc=%d segments=%zu\n",rc,c.segment_count);
    check(rc==GEOMLIB_MESH_SUCCESS,"contour extraction succeeds");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(c.segment_count>0,"iso line has segments");
        double err=0.0;
        bool vertical=true;
        for(size_t i=0;i<c.segment_count;++i){
            geomlib_vec3 a=c.segments[i*2],b=c.segments[i*2+1];
            err=std::fmax(err,std::fabs(a.x-2.0));
            err=std::fmax(err,std::fabs(b.x-2.0));
            if(std::fabs(a.x-b.x)>1e-4) vertical=false;
        }
        std::printf("contour max value error=%.6f vertical=%d\n",err,vertical);
        check(err<1e-5,"all endpoints lie on the iso value");
        check(vertical,"isolines are vertical for the x field");
        geomlib_mesh_contour_free(&c);
    }

    rc=geomlib_mesh_contour_lines(&m,vals.data(),100.0f,&c);
    check(rc==GEOMLIB_MESH_SUCCESS && c.segment_count==0,
          "out-of-range iso gives no segments");
    geomlib_mesh_contour_free(&c);
    geomlib_trimesh_free(&m);
    if(failures!=0){std::printf("%d contour test(s) failed\n",failures);return 1;}
    std::printf("All contour tests passed\n");
    return 0;
}
