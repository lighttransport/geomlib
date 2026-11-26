/*
 * Rasterization: a square's inside count approximates its area and signed
 * distances are negative inside.
 */
#include <geomlib/geomlib_raster2d.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_vec2 sq[4]={
        geomlib_vec2_create(0,0),geomlib_vec2_create(2,0),
        geomlib_vec2_create(2,2),geomlib_vec2_create(0,2)
    };
    geomlib_raster2d r=geomlib_raster2d_create();
    int rc=geomlib_rasterize_polygon2d(sq,4,50,50,0.1f,1,&r);
    size_t inside=0;
    for(size_t i=0;i<(size_t)r.nx*r.ny;++i) if(r.mask[i]) inside++;
    double area=(double)inside*r.cell_size*r.cell_size;
    std::printf("raster rc=%d inside=%zu area=%.4f expect=4\n",rc,inside,area);
    check(rc==GEOMLIB_MESH_SUCCESS,"rasterization succeeds");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(area-4.0)<0.1,"area approximates square");
        bool signs=true;
        for(size_t i=0;i<(size_t)r.nx*r.ny;++i) {
            if(r.mask[i] && r.distance[i]>0) signs=false;
            if(!r.mask[i] && r.distance[i]<0) signs=false;
        }
        check(signs,"signed distances match mask");
        geomlib_raster2d_free(&r);
    }
    if(failures!=0){std::printf("%d raster test(s) failed\n",failures);return 1;}
    std::printf("All raster tests passed\n");
    return 0;
}
