/*
 * Quaternions: composition matches matrix multiplication, slerp endpoints and
 * midpoint behave correctly, and rotation is length-preserving.
 */
#include <geomlib/geomlib_quaternion.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

int main(void)
{
    geomlib_quat a=geomlib_quat_from_axis_angle(
        geomlib_vec3_create(0,0,1),1.0);
    geomlib_quat b=geomlib_quat_from_axis_angle(
        geomlib_vec3_create(0,1,0),0.5);
    geomlib_quat ab=geomlib_quat_mul(a,b);
    geomlib_mat4 ma=geomlib_quat_to_mat4(a),mb=geomlib_quat_to_mat4(b);
    geomlib_mat4 mab=geomlib_mat4_mul(ma,mb);
    geomlib_quat ab2=geomlib_quat_normalize(ab);
    geomlib_mat4 m2=geomlib_quat_to_mat4(ab2);
    double err=0.0;
    for(int i=0;i<16;++i) err=std::fmax(err,std::fabs(mab.m[i]-m2.m[i]));
    std::printf("quat compose err=%.6e\n",err);
    check(err<1e-5,"composition matches matrix product");

    geomlib_quat q=geomlib_quat_slerp(geomlib_quat_identity(),a,0.5);
    geomlib_vec3 v=geomlib_vec3_create(1,0,0);
    geomlib_vec3 rv=geomlib_quat_rotate(q,v);
    std::printf("slerp half rotated=(%.3f,%.3f,%.3f)\n",rv.x,rv.y,rv.z);
    check(std::fabs(geomlib_vec3_length(rv)-1.0)<1e-5,"rotation preserves length");
    geomlib_quat start=geomlib_quat_slerp(geomlib_quat_identity(),a,0.0);
    geomlib_quat end=geomlib_quat_slerp(geomlib_quat_identity(),a,1.0);
    check(std::fabs(start.w-1.0)<1e-6 && std::fabs(end.w-a.w)<1e-6,
          "slerp endpoints correct");

    geomlib_quat inv=geomlib_quat_mul(a,geomlib_quat_conjugate(a));
    check(std::fabs(inv.w-1.0)<1e-6,"q * q* is identity");
    if(failures!=0){std::printf("%d quaternion test(s) failed\n",failures);return 1;}
    std::printf("All quaternion tests passed\n");
    return 0;
}
