/*
 * Horn rigid fitting: recovers a known rotation + translation, and a known
 * uniform scale when requested.
 */
#include <geomlib/geomlib_rigid_fit.h>

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
    std::vector<geomlib_vec3> src,dst;
    for(int i=0;i<40;++i){
        float x=(float)((i*37)%101)/10.0f;
        float y=(float)((i*53)%97)/10.0f;
        float z=(float)((i*29)%89)/10.0f;
        src.push_back(geomlib_vec3_create(x,y,z));
    }
    geomlib_mat4 known=geomlib_mat4_rotate_axis(
        geomlib_vec3_create(0.3f,0.5f,0.8f),0.7f);
    known=geomlib_mat4_mul(geomlib_mat4_translate(1.5f,-2.0f,0.75f),known);
    for(auto &p:src){
        geomlib_vec4 v=geomlib_mat4_mul_vec4(known,
            geomlib_vec4_create(p.x,p.y,p.z,1.0f));
        dst.push_back(geomlib_vec3_create(v.x,v.y,v.z));
    }
    geomlib_mat4 fit;
    float rms=-1.0f;
    int rc=geomlib_fit_rigid_points(src.data(),dst.data(),src.size(),false,&fit,&rms);
    std::printf("rigid fit rc=%d rms=%.6f\n",rc,rms);
    check(rc==GEOMLIB_MESH_SUCCESS,"rigid fit succeeds");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(rms<1e-3f,"rigid fit recovers transform");
        double err=0.0;
        for(size_t i=0;i<src.size();++i){
            geomlib_vec4 v=geomlib_mat4_mul_vec4(fit,
                geomlib_vec4_create(src[i].x,src[i].y,src[i].z,1.0f));
            err+=std::hypot((double)v.x-dst[i].x,(double)v.y-dst[i].y);
        }
        check(err<1e-2,"fit transform reproduces points");
    }

    /* Uniform scale 2 plus the same rotation/translation. */
    geomlib_mat4 scaled=geomlib_mat4_mul(known,geomlib_mat4_scale_xyz(2.0f,2.0f,2.0f));
    for(size_t i=0;i<src.size();++i){
        geomlib_vec4 v=geomlib_mat4_mul_vec4(scaled,
            geomlib_vec4_create(src[i].x,src[i].y,src[i].z,1.0f));
        dst[i]=geomlib_vec3_create(v.x,v.y,v.z);
    }
    rc=geomlib_fit_rigid_points(src.data(),dst.data(),src.size(),true,&fit,&rms);
    std::printf("scale fit rc=%d rms=%.6f\n",rc,rms);
    check(rc==GEOMLIB_MESH_SUCCESS && rms<1e-3f,
          "similarity fit recovers uniform scale");
    if(failures!=0){std::printf("%d rigid fit test(s) failed\n",failures);return 1;}
    std::printf("All rigid fit tests passed\n");
    return 0;
}
