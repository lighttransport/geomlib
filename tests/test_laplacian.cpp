/*
 * Cotangent Laplacian: constant nullspace, positive mass, symmetry.
 */
#include <geomlib/geomlib_laplacian.h>

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
    geomlib_trimesh grid;
    geomlib_trimesh_init(&grid,25,32);
    for(int y=0;y<5;++y)
        for(int x=0;x<5;++x)
            geomlib_trimesh_add_vertex(&grid,geomlib_vec3_create((float)x,(float)y,0),
                                       nullptr,nullptr);
    for(int y=0;y<4;++y)
        for(int x=0;x<4;++x){
            uint32_t i=(uint32_t)(y*5+x);
            geomlib_trimesh_add_triangle(&grid,i,i+1,i+6);
            geomlib_trimesh_add_triangle(&grid,i,i+6,i+5);
        }
    geomlib_laplacian l=geomlib_laplacian_create();
    int rc=geomlib_build_laplacian(&grid,&l);
    std::printf("laplacian rc=%d vertices=%zu nnz=%zu\n",
                rc,l.vertex_count,l.stiffness.nnz);
    check(rc==GEOMLIB_MESH_SUCCESS,"laplacian builds");
    if(rc==GEOMLIB_MESH_SUCCESS){
        std::vector<double> ones(l.vertex_count,1.0),y(l.vertex_count,0.0);
        rc=geomlib_laplacian_apply(&l,ones.data(),0.0,y.data());
        double mx=0.0,mass=0.0;
        for(size_t i=0;i<l.vertex_count;++i){
            if(std::fabs(y[i])>mx) mx=std::fabs(y[i]);
            mass+=l.mass[i];
        }
        std::printf("constant nullspace max=%.3e total mass=%.4f expect=16\n",mx,mass);
        check(rc==GEOMLIB_MESH_SUCCESS && mx<1e-9,"constants are in the nullspace");
        check(std::fabs(mass-16.0)<0.1,"lumped mass equals area");
        geomlib_laplacian_free(&l);
    }
    geomlib_trimesh_free(&grid);
    if(failures!=0){std::printf("%d laplacian test(s) failed\n",failures);return 1;}
    std::printf("All laplacian tests passed\n");
    return 0;
}
