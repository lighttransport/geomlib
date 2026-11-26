/*
 * Point-to-plane ICP: realign a perturbed sphere onto its original position.
 */
#include <geomlib/geomlib_icp_plane.h>
#include <geomlib/geomlib_bvh.h>

#include <cmath>
#include <cstdio>

static int failures=0;
static void check(bool ok,const char*name)
{
    if(!ok){std::printf("FAIL %s\n",name);failures++;}
}

static float rms_to_target(const geomlib_trimesh &s, const geomlib_bvh &bvh)
{
    double sum=0.0;
    for(size_t i=0;i<s.vertex_count;++i){
        geomlib_vec3 cp; float d=0.0f;
        if(!geomlib_bvh_closest_point(&bvh,s.vertices[i],&cp,&d,nullptr))
            d=10.0f;
        sum+=d*d;
    }
    return (float)std::sqrt(sum/s.vertex_count);
}

int main(void)
{
    geomlib_trimesh target={},source={};
    int rc=geomlib_create_sphere(&target,1.0f,16,12);
    check(rc==GEOMLIB_MESH_SUCCESS,"target sphere");
    if(rc!=GEOMLIB_MESH_SUCCESS) return 1;
    geomlib_trimesh_init(&source,target.vertex_count,target.triangle_count);
    for(size_t i=0;i<target.vertex_count;++i)
        geomlib_trimesh_add_vertex(&source,
            geomlib_vec3_create(target.vertices[i].x+0.15f,
                                target.vertices[i].y-0.10f,
                                target.vertices[i].z+0.05f),nullptr,nullptr);
    for(size_t i=0;i<target.triangle_count;++i)
        geomlib_trimesh_add_triangle(&source,target.indices[i*3],
            target.indices[i*3+1],target.indices[i*3+2]);

    geomlib_bvh bvh;
    geomlib_bvh_build(&bvh,&target,GEOMLIB_BVH_QUALITY_MEDIUM);
    float before=rms_to_target(source,bvh);

    geomlib_mat4 T=geomlib_mat4_identity();
    float rms=0.0f;
    rc=geomlib_icp_point_to_plane(&source,&target,nullptr,&T,&rms);
    std::printf("icp plane rc=%d rms=%.5f before=%.5f\n",rc,rms,before);
    check(rc==GEOMLIB_MESH_SUCCESS,"ICP succeeds");

    /* Apply the fitted transform and re-measure. */
    for(size_t i=0;i<source.vertex_count;++i){
        geomlib_vec4 v=geomlib_mat4_mul_vec4(T,
            geomlib_vec4_create(source.vertices[i].x,source.vertices[i].y,
                                source.vertices[i].z,1.0f));
        source.vertices[i]=geomlib_vec3_create(v.x,v.y,v.z);
    }
    float after=rms_to_target(source,bvh);
    std::printf("icp plane after=%.5f\n",after);
    check(after<before*0.2f,"point-to-plane ICP reduces misalignment");
    check(std::fabs(after-rms)<1e-3f,"reported RMS matches measured RMS");

    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&source);
    geomlib_trimesh_free(&target);
    if(failures!=0){std::printf("%d ICP plane test(s) failed\n",failures);return 1;}
    std::printf("All ICP plane tests passed\n");
    return 0;
}
