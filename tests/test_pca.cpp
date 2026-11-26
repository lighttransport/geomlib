/*
 * PCA: an anisotropic point cloud's principal axis follows its long direction
 * and eigenvalues are ordered descending.
 */
#include <geomlib/geomlib_pca.h>

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
    std::vector<geomlib_vec3> pts;
    for(int i=0;i<200;++i){
        float u=(float)((i*37)%101)/10.0f-5.0f;
        float v=(float)((i*53)%97)/100.0f-0.5f;
        pts.push_back(geomlib_vec3_create(u*3.0f, v, 0.1f*v));
    }
    geomlib_pca_result p;
    int rc=geomlib_compute_pca(pts.data(),pts.size(),&p);
    std::printf("pca rc=%d mean=(%.3f,%.3f,%.3f) eig=(%.3f,%.3f,%.3f)\n",
                rc,p.mean.x,p.mean.y,p.mean.z,
                p.eigenvalues[0],p.eigenvalues[1],p.eigenvalues[2]);
    check(rc==GEOMLIB_MESH_SUCCESS,"PCA computes");
    if(rc==GEOMLIB_MESH_SUCCESS){
        check(p.eigenvalues[0]>=p.eigenvalues[1] &&
              p.eigenvalues[1]>=p.eigenvalues[2],"eigenvalues descending");
        check(std::fabs(p.eigenvectors[0].x)>0.99,"first axis is x");
        check(std::fabs(p.eigenvalues[0])>10.0*p.eigenvalues[1],
              "long direction dominates");
    }
    if(failures!=0){std::printf("%d PCA test(s) failed\n",failures);return 1;}
    std::printf("All PCA tests passed\n");
    return 0;
}
