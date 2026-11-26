/*
 * Harmonic cage: weights sum to one, identity is a no-op and an affine cage
 * deformation is reproduced exactly.
 */
#include <geomlib/geomlib_cage2d_harmonic.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_vec2 cage[4]={
        geomlib_vec2_create(0,0),geomlib_vec2_create(2,0),
        geomlib_vec2_create(2,2),geomlib_vec2_create(0,2)
    };
    geomlib_vec2 q=geomlib_vec2_create(0.7f,0.4f);
    double w[4];
    int rc=geomlib_cage2d_harmonic_weights(cage,4,q,w);
    double sum=0.0;
    for(int i=0;i<4;++i) sum+=w[i];
    std::printf("cage harmonic rc=%d sum=%.6f w=(%.3f,%.3f,%.3f,%.3f)\n",
                rc,sum,w[0],w[1],w[2],w[3]);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(sum-1.0)<1e-6,
          "weights sum to one");

    geomlib_vec2 out;
    rc=geomlib_cage2d_harmonic_deform(cage,cage,4,&q,1,&out);
    check(rc==GEOMLIB_MESH_SUCCESS &&
          geomlib_vec2_length(geomlib_vec2_sub(out,q))<1e-4f,
          "identity cage leaves point unchanged");

    geomlib_vec2 moved[4];
    for(int i=0;i<4;++i)
        moved[i]=geomlib_vec2_create(cage[i].x*2.0f+1.0f,cage[i].y*0.5f-1.0f);
    rc=geomlib_cage2d_harmonic_deform(cage,moved,4,&q,1,&out);
    geomlib_vec2 expect=geomlib_vec2_create(
        q.x*2.0f+1.0f,q.y*0.5f-1.0f);
    std::printf("affine err=%.6f\n",geomlib_vec2_length(geomlib_vec2_sub(out,expect)));
    check(rc==GEOMLIB_MESH_SUCCESS &&
          geomlib_vec2_length(geomlib_vec2_sub(out,expect))<1e-3f,
          "affine cage transform reproduced");
    if(failures!=0){std::printf("%d harmonic cage test(s) failed\n",failures);return 1;}
    std::printf("All harmonic cage tests passed\n");
    return 0;
}
