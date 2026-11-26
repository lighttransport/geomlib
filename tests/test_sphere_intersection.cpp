/*
 * Sphere-sphere intersection: analytic circle, tangency, disjoint and
 * containment cases.
 */
#include <geomlib/geomlib_sphere_intersection.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_sphere_intersection r=geomlib_sphere_sphere_intersection(
        geomlib_vec3_create(0,0,0),2.0f,geomlib_vec3_create(3,0,0),2.0f);
    std::printf("sphere int kind=%d center=%.3f radius=%.4f expect=1.323\n",
                (int)r.kind,r.center.x,r.radius);
    check(r.kind==GEOMLIB_SPHERE_CIRCLE && std::fabs(r.center.x-1.5)<1e-4 &&
          std::fabs(r.radius-std::sqrt(4.0-2.25))<1e-3,
          "unit-offset intersection circle is analytic");

    r=geomlib_sphere_sphere_intersection(
        geomlib_vec3_create(0,0,0),1.0f,geomlib_vec3_create(2,0,0),1.0f);
    check(r.kind==GEOMLIB_SPHERE_POINT && std::fabs(r.center.x-1.0)<1e-4,
          "externally tangent at one point");

    r=geomlib_sphere_sphere_intersection(
        geomlib_vec3_create(0,0,0),1.0f,geomlib_vec3_create(5,0,0),1.0f);
    check(r.kind==GEOMLIB_SPHERE_NONE,"disjoint spheres");

    r=geomlib_sphere_sphere_intersection(
        geomlib_vec3_create(0,0,0),1.0f,geomlib_vec3_create(0,0,0),1.0f);
    check(r.kind==GEOMLIB_SPHERE_CONTAINED,"coincident spheres");
    if(failures!=0){std::printf("%d sphere intersection test(s) failed\n",failures);return 1;}
    std::printf("All sphere intersection tests passed\n");
    return 0;
}
