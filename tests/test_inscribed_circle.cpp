/*
 * Largest inscribed circle: square and rectangle give the analytic inradius,
 * and the circle touches at least two edges.
 */
#include <geomlib/geomlib_inscribed_circle.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static double seg_dist(geomlib_vec2 q, geomlib_vec2 a, geomlib_vec2 b)
{
    double vx=b.x-a.x,vy=b.y-a.y,wx=q.x-a.x,wy=q.y-a.y;
    double l2=vx*vx+vy*vy,t=l2>0?(wx*vx+wy*vy)/l2:0;
    if(t<0)t=0; if(t>1)t=1;
    return std::hypot((double)q.x-(a.x+t*(b.x-a.x)),
                      (double)q.y-(a.y+t*(b.y-a.y)));
}

int main(void)
{
    geomlib_vec2 sq[4]={
        geomlib_vec2_create(-1,-1),geomlib_vec2_create(1,-1),
        geomlib_vec2_create(1,1),geomlib_vec2_create(-1,1)
    };
    geomlib_inscribed_circle c;
    int rc=geomlib_largest_inscribed_circle(sq,4,0.05f,&c);
    std::printf("circle square rc=%d center=(%.3f,%.3f) radius=%.4f expect=1\n",
                rc,c.center.x,c.center.y,c.radius);
    check(rc==GEOMLIB_MESH_SUCCESS,"square circle computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(c.radius-1.0f)<0.02f,"square inradius is 1");
        check(std::hypot(c.center.x,c.center.y)<0.02f,"square centre is origin");
        double d[4];
        for(int k=0;k<4;++k) d[k]=seg_dist(c.center,sq[k],sq[(k+1)%4]);
        int touches=0;
        for(int k=0;k<4;++k) if(std::fabs(d[k]-c.radius)<0.02) touches++;
        check(touches>=2,"circle touches at least two edges");
    }

    geomlib_vec2 rect[4]={
        geomlib_vec2_create(-2,-1),geomlib_vec2_create(2,-1),
        geomlib_vec2_create(2,1),geomlib_vec2_create(-2,1)
    };
    rc=geomlib_largest_inscribed_circle(rect,4,0.05f,&c);
    std::printf("circle rect rc=%d radius=%.4f expect=1\n",rc,c.radius);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(c.radius-1.0f)<0.02f,
          "rectangle inradius is 1");
    if(failures!=0){std::printf("%d inscribed circle test(s) failed\n",failures);return 1;}
    std::printf("All inscribed circle tests passed\n");
    return 0;
}
