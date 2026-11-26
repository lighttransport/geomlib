/*
 * Tetrahedral cage: barycentric weights sum to 1, identity cage leaves
 * points unchanged, and an affine cage transform is reproduced exactly.
 */
#include <geomlib/geomlib_cage3d.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_vec3 cage[4]={
        geomlib_vec3_create(0,0,0),geomlib_vec3_create(1,0,0),
        geomlib_vec3_create(0,1,0),geomlib_vec3_create(0,0,1)
    };
    geomlib_vec3 q=geomlib_vec3_create(0.2f,0.3f,0.1f);
    double w[4];
    int rc=geomlib_tetra_weights(cage,q,w);
    double sum=0.0;
    for(int i=0;i<4;++i) sum+=w[i];
    std::printf("cage weights rc=%d sum=%.6f w=(%.3f,%.3f,%.3f,%.3f)\n",
                rc,sum,w[0],w[1],w[2],w[3]);
    check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(sum-1.0)<1e-6,
          "weights sum to one");

    geomlib_vec3 out;
    rc=geomlib_deform_tetra_cage(cage,cage,&q,1,&out);
    check(rc==GEOMLIB_MESH_SUCCESS &&
          geomlib_vec3_length(geomlib_vec3_sub(out,q))<1e-5f,
          "identity cage leaves points unchanged");

    /* Affine cage deformation: scale by 2 and translate by (1,2,3). */
    geomlib_vec3 moved[4];
    for(int i=0;i<4;++i)
        moved[i]=geomlib_vec3_add(geomlib_vec3_mul(cage[i],2.0f),
                                  geomlib_vec3_create(1,2,3));
    rc=geomlib_deform_tetra_cage(cage,moved,&q,1,&out);
    geomlib_vec3 expect=geomlib_vec3_add(geomlib_vec3_mul(q,2.0f),
                                          geomlib_vec3_create(1,2,3));
    std::printf("cage affine err=%.6f\n",
                geomlib_vec3_length(geomlib_vec3_sub(out,expect)));
    check(rc==GEOMLIB_MESH_SUCCESS &&
          geomlib_vec3_length(geomlib_vec3_sub(out,expect))<1e-4f,
          "affine cage transform is reproduced exactly");
    if(failures!=0){std::printf("%d cage test(s) failed\n",failures);return 1;}
    std::printf("All cage tests passed\n");
    return 0;
}
