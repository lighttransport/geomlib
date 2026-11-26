/*
 * Boolean operations on polygons with holes: a framed square unioned with a
 * small square inside its hole gives the frame plus the square.
 */
#include <geomlib/geomlib_polybool_holes.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static double area(const geomlib_poly_bool_result &r)
{
    double s=0.0;
    for(size_t l=0;l<r.loop_count;++l){
        size_t a=r.offsets[l],b=r.offsets[l+1];
        double t=0.0;
        for(size_t i=a;i<b;++i){
            size_t j=(i+1<b)?i+1:a;
            t+=(double)r.loops[i].x*r.loops[j].y
             -(double)r.loops[j].x*r.loops[i].y;
        }
        s+=std::fabs(0.5*t);
    }
    return s;
}

int main(void)
{
    geomlib_vec2 outer[4]={
        geomlib_vec2_create(0,0),geomlib_vec2_create(4,0),
        geomlib_vec2_create(4,4),geomlib_vec2_create(0,4)
    };
    geomlib_vec2 hole[4]={
        geomlib_vec2_create(1,1),geomlib_vec2_create(3,1),
        geomlib_vec2_create(3,3),geomlib_vec2_create(1,3)
    };
    const geomlib_vec2 *holes[1]={hole};
    size_t sizes[1]={4};
    geomlib_polygon_with_holes A;
    int rc=geomlib_polygon_with_holes_set(outer,4,holes,sizes,1,&A);
    check(rc==GEOMLIB_MESH_SUCCESS,"framed polygon built");

    geomlib_polygon_with_holes B;
    geomlib_vec2 inner[4]={
        geomlib_vec2_create(1.5f,1.5f),geomlib_vec2_create(2.5f,1.5f),
        geomlib_vec2_create(2.5f,2.5f),geomlib_vec2_create(1.5f,2.5f)
    };
    rc=geomlib_polygon_with_holes_set(inner,4,nullptr,nullptr,0,&B);
    check(rc==GEOMLIB_MESH_SUCCESS,"inner square built");

    geomlib_poly_bool_result r=geomlib_poly_bool_result_create();
    rc=geomlib_polygon_boolean_holes(&A,&B,GEOMLIB_POLY_UNION,&r);
    std::printf("union rc=%d loops=%zu area=%.4f expect=13\n",rc,r.loop_count,area(r));
    check(rc==GEOMLIB_MESH_SUCCESS,"union computes");
    if(rc==GEOMLIB_MESH_SUCCESS)
        check(std::fabs(area(r)-13.0)<1e-3,"union area is exactly 13");
    geomlib_poly_bool_result_free(&r);

    r=geomlib_poly_bool_result_create();
    rc=geomlib_polygon_boolean_holes(&A,&B,GEOMLIB_POLY_INTERSECTION,&r);
    std::printf("intersection rc=%d loops=%zu area=%.4f expect=0\n",
                rc,r.loop_count,area(r));
    check(rc==GEOMLIB_MESH_SUCCESS && area(r)<0.1,
          "intersection with hole is empty");
    geomlib_poly_bool_result_free(&r);

    geomlib_polygon_with_holes_free(&A);
    geomlib_polygon_with_holes_free(&B);
    if(failures!=0){std::printf("%d polybool holes test(s) failed\n",failures);return 1;}
    std::printf("All polybool holes tests passed\n");
    return 0;
}
