/*
 * Mass properties: a unit box has volume 1, centre (0.5,0.5,0.5) and inertia
 * 1/6 about its centre; a sphere has the analytic inertia 8*pi/15.
 */
#include <geomlib/geomlib_mass_properties.h>

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
    int rc=geomlib_create_box(&box,1.0f,1.0f,1.0f);
    check(rc==GEOMLIB_MESH_SUCCESS,"box created");
    if(rc!=GEOMLIB_MESH_SUCCESS) return 1;
    geomlib_mass_properties mp;
    rc=geomlib_compute_mass_properties(&box,&mp);
    std::printf("mass box rc=%d vol=%.4f center=(%.3f,%.3f,%.3f) I=(%.4f,%.4f,%.4f)\n",
                rc,mp.volume,mp.center.x,mp.center.y,mp.center.z,
                mp.inertia.m[0],mp.inertia.m[4],mp.inertia.m[8]);
    check(rc==GEOMLIB_MESH_SUCCESS,"mass properties compute");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(std::fabs(mp.volume-1.0)<1e-3,"box volume is 1");
        check(std::fabs(mp.center.x)<1e-3 && std::fabs(mp.center.y)<1e-3 &&
              std::fabs(mp.center.z)<1e-3,"box centre is at the origin");
        check(std::fabs(mp.inertia.m[0]-1.0/6.0)<0.01 &&
              std::fabs(mp.inertia.m[4]-1.0/6.0)<0.01 &&
              std::fabs(mp.inertia.m[8]-1.0/6.0)<0.01,
              "box inertia is 1/6");
    }
    geomlib_trimesh_free(&box);

    geomlib_trimesh sphere={};
    rc=geomlib_create_sphere(&sphere,1.0f,24,16);
    if(rc==GEOMLIB_MESH_SUCCESS){
        rc=geomlib_compute_mass_properties(&sphere,&mp);
        std::printf("mass sphere rc=%d vol=%.4f expect=%.4f I=%.4f expect=%.4f\n",
                    rc,mp.volume,4*M_PI/3,mp.inertia.m[0],8*M_PI/15);
        check(rc==GEOMLIB_MESH_SUCCESS && std::fabs(mp.volume-4*M_PI/3)<0.1,
              "sphere volume is analytic");
        if(rc==GEOMLIB_MESH_SUCCESS)
            check(std::fabs(mp.inertia.m[0]-8*M_PI/15)<0.15,
                  "sphere inertia is analytic");
        geomlib_trimesh_free(&sphere);
    }
    if(failures!=0){std::printf("%d mass property test(s) failed\n",failures);return 1;}
    std::printf("All mass property tests passed\n");
    return 0;
}
