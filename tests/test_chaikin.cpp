/*
 * Chaikin subdivision: endpoints preserved, point count doubles per level,
 * and a straight line stays straight.
 */
#include <geomlib/geomlib_chaikin.h>

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
    std::vector<geomlib_vec3> line={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(3,0,0),
        geomlib_vec3_create(6,0,0)
    };
    std::vector<geomlib_vec3> out(9);
    size_t n=geomlib_chaikin_subdivide(line.data(),line.size(),2,out.data(),out.size());
    std::printf("chaikin line rc=%zu first_x=%.3f last_x=%.3f\n",
                n,out[0].x,out[n-1].x);
    check(n==9,"two levels give 9 points");
    if(n==9){
        bool straight=true;
        for(size_t i=0;i<n;++i)
            if(std::fabs(out[i].y)>1e-6f||std::fabs(out[i].z)>1e-6f) straight=false;
        check(straight,"straight line stays straight");
        check(std::fabs(out[0].x)<1e-6f && std::fabs(out[n-1].x-6.0f)<1e-6f,
              "endpoints preserved");
    }

    std::vector<geomlib_vec3> corner={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,1,0),
        geomlib_vec3_create(2,0,0)
    };
    std::vector<geomlib_vec3> refined(5);
    n=geomlib_chaikin_subdivide(corner.data(),corner.size(),1,refined.data(),refined.size());
    std::printf("chaikin corner rc=%zu mid=(%.3f,%.3f)\n",n,refined[2].x,refined[2].y);
    check(n==5,"one level gives 5 points");
    check(refined[0].x==0.0f && refined[4].x==2.0f,"corner endpoints preserved");
    check(refined[2].y>0.5f,"corner apex is cut but retained");
    if(failures!=0){std::printf("%d Chaikin test(s) failed\n",failures);return 1;}
    std::printf("All Chaikin tests passed\n");
    return 0;
}
