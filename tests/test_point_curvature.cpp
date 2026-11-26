/*
 * Local PCA curvature: planar points have near-zero curvature and vertical
 * normals; sphere points have higher curvature and radial normals.
 */
#include <geomlib/geomlib_point_curvature.h>

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
    geomlib_pointcloud plane;
    geomlib_pointcloud_init(&plane,300);
    for(int i=0;i<300;++i)
        geomlib_pointcloud_add_point(&plane,
            geomlib_vec3_create(((i*37)%101)/100.0f,((i*73)%97)/100.0f,0),
            nullptr,1.0f);
    std::vector<float> pc(plane.point_count);
    std::vector<geomlib_vec3> pn(plane.point_count);
    int rc=geomlib_pointcloud_curvature(&plane,20,pc.data(),pn.data());
    float mx=0.0f;
    for(auto c:pc) if(c>mx) mx=c;
    std::printf("curvature plane rc=%d max=%.6f\n",rc,mx);
    check(rc==GEOMLIB_MESH_SUCCESS,"plane curvature succeeds");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(mx<0.01f,"planar curvature is near zero");
        bool vertical=true;
        for(auto &n:pn) if(std::fabs(n.z)<0.95f) vertical=false;
        check(vertical,"planar normals are vertical");
    }
    geomlib_pointcloud_free(&plane);

    geomlib_pointcloud sphere;
    geomlib_pointcloud_init(&sphere,400);
    for(int i=0;i<400;++i){
        float u=i*0.618033f*6.283185f;
        float z=1.0f-2.0f*(float)(i+0.5)/400.0f;
        float r=std::sqrt(std::fmax(0.0f,1.0f-z*z));
        geomlib_pointcloud_add_point(&sphere,
            geomlib_vec3_create(r*std::cos(u),r*std::sin(u),z),nullptr,1.0f);
    }
    std::vector<float> sc(sphere.point_count);
    std::vector<geomlib_vec3> sn(sphere.point_count);
    rc=geomlib_pointcloud_curvature(&sphere,20,sc.data(),sn.data());
    float avg=0.0f;
    for(auto c:sc) avg+=c;
    avg/=sc.size();
    float radial=0.0f;
    for(size_t i=0;i<sphere.point_count;++i)
        radial+=std::fabs(geomlib_vec3_dot(sn[i],geomlib_vec3_normalize(sphere.points[i])));
    radial/=(float)sphere.point_count;
    std::printf("curvature sphere rc=%d avg=%.5f radial_dot=%.3f\n",rc,avg,radial);
    check(rc==GEOMLIB_MESH_SUCCESS,"sphere curvature succeeds");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(avg>0.01f,"sphere curvature is positive");
        check(radial>0.8f,"sphere normals point radially");
    }
    geomlib_pointcloud_free(&sphere);
    if(failures!=0){std::printf("%d point curvature test(s) failed\n",failures);return 1;}
    std::printf("All point curvature tests passed\n");
    return 0;
}
