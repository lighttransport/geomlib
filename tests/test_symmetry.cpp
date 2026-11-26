/*
 * Symmetry detection: a box is symmetric (near-zero RMS), an asymmetric
 * point-cloud-like mesh has a larger error.
 */
#include <geomlib/geomlib_symmetry.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_trimesh box={};
    int rc=geomlib_create_box(&box,2.0f,1.0f,0.5f);
    check(rc==GEOMLIB_MESH_SUCCESS,"box created");
    if(rc!=GEOMLIB_MESH_SUCCESS) return 1;
    geomlib_symmetry_result s;
    rc=geomlib_detect_symmetry(&box,&s);
    std::printf("symmetry box rc=%d rms=%.5f normal=(%.2f,%.2f,%.2f)\n",
                rc,s.rms_error,s.plane_normal.x,s.plane_normal.y,s.plane_normal.z);
    check(rc==GEOMLIB_MESH_SUCCESS,"symmetry computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(s.rms_error<0.02f,"box is symmetric");
        check(std::fabs(geomlib_vec3_length(s.plane_normal)-1.0f)<1e-3f,
              "plane normal is unit length");
    }
    geomlib_trimesh_free(&box);

    /* Perturb one vertex far away: symmetry error grows. */
    geomlib_trimesh asym={};
    rc=geomlib_create_box(&asym,2.0f,1.0f,0.5f);
    asym.vertices[1].x+=1.5f;
    rc=geomlib_detect_symmetry(&asym,&s);
    std::printf("symmetry asymmetric rc=%d rms=%.5f\n",rc,s.rms_error);
    check(rc==GEOMLIB_MESH_SUCCESS && s.rms_error>0.02f,
          "perturbed box is less symmetric");
    geomlib_trimesh_free(&asym);
    if(failures!=0){std::printf("%d symmetry test(s) failed\n",failures);return 1;}
    std::printf("All symmetry tests passed\n");
    return 0;
}
