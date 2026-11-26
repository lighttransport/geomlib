/*
 * Convex decomposition: an L-shape splits into convex pieces whose total
 * area equals the original polygon.
 */
#include <geomlib/geomlib_convex_decompose2d.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static double cross(const geomlib_vec2 &a,const geomlib_vec2 &b,const geomlib_vec2 &c)
{
    return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
}

static bool convex(const geomlib_vec2 *p,size_t n)
{
    bool sign=true;
    for(size_t i=0;i<n;++i){
        double c=cross(p[i],p[(i+1)%n],p[(i+2)%n]);
        if(c<-1e-8) sign=false;
    }
    /* Require consistent (all convex) winding: no strictly reflex turn. */
    bool neg=true;
    for(size_t i=0;i<n;++i){
        double c=cross(p[i],p[(i+1)%n],p[(i+2)%n]);
        if(c>1e-8) neg=false;
    }
    return sign||neg;
}

static double area(const geomlib_vec2 *p,size_t n)
{
    double s=0.0;
    for(size_t i=0;i<n;++i){
        size_t j=(i+1)%n;
        s+=p[i].x*p[j].y-p[j].x*p[i].y;
    }
    return 0.5*std::fabs(s);
}

int main(void)
{
    geomlib_vec2 L[6]={
        geomlib_vec2_create(0,0),geomlib_vec2_create(2,0),
        geomlib_vec2_create(2,1),geomlib_vec2_create(1,1),
        geomlib_vec2_create(1,2),geomlib_vec2_create(0,2)
    };
    geomlib_poly_decomposition d=geomlib_poly_decomposition_create();
    int rc=geomlib_convex_decompose2d(L,6,&d);
    double total=0.0;
    bool all_convex=true;
    for(size_t k=0;k<d.part_count;++k){
        size_t a=d.offsets[k],b=d.offsets[k+1];
        total+=area(&d.parts[a],b-a);
        if(!convex(&d.parts[a],b-a)) all_convex=false;
    }
    std::printf("decompose rc=%d parts=%zu area=%.4f expect=3\n",
                rc,d.part_count,total);
    check(rc==GEOMLIB_MESH_SUCCESS,"decomposition succeeds");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(d.part_count>=2,"L shape splits into multiple parts");
        check(all_convex,"all parts are convex");
        check(std::fabs(total-3.0)<1e-4,"total area preserved");
    }
    geomlib_poly_decomposition_free(&d);

    geomlib_vec2 sq[4]={
        geomlib_vec2_create(0,0),geomlib_vec2_create(1,0),
        geomlib_vec2_create(1,1),geomlib_vec2_create(0,1)
    };
    rc=geomlib_convex_decompose2d(sq,4,&d);
    std::printf("decompose square rc=%d parts=%zu\n",rc,d.part_count);
    check(rc==GEOMLIB_MESH_SUCCESS && d.part_count==1,"square stays one part");
    geomlib_poly_decomposition_free(&d);

    if(failures!=0){std::printf("%d decompose test(s) failed\n",failures);return 1;}
    std::printf("All decompose tests passed\n");
    return 0;
}
