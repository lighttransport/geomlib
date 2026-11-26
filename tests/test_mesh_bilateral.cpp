/*
 * Bilateral mesh smoothing: reduces height noise on a planar patch while
 * keeping vertex normals unit and all positions finite.
 */
#include <geomlib/geomlib_mesh_bilateral.h>

#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static geomlib_trimesh make_grid()
{
    geomlib_trimesh m;
    geomlib_trimesh_init(&m,25,32);
    for (int y=0;y<5;++y)
        for (int x=0;x<5;++x)
            geomlib_trimesh_add_vertex(&m,geomlib_vec3_create((float)x,(float)y,0),
                                       nullptr,nullptr);
    for (int y=0;y<4;++y)
        for (int x=0;x<4;++x) {
            uint32_t i=(uint32_t)(y*5+x);
            geomlib_trimesh_add_triangle(&m,i,i+1,i+6);
            geomlib_trimesh_add_triangle(&m,i,i+6,i+5);
        }
    return m;
}

int main(void)
{
    geomlib_trimesh flat=make_grid();
    int rc=geomlib_mesh_bilateral_smooth(&flat,nullptr);
    float maxz=0.0f;
    for (size_t i=0;i<flat.vertex_count;++i) maxz=std::fmax(maxz,std::fabs(flat.vertices[i].z));
    std::printf("bilateral flat rc=%d max |z|=%.6f\n", rc, maxz);
    check(rc==GEOMLIB_MESH_SUCCESS,"flat smoothing succeeds");
    check(maxz<1e-5f,"planar mesh stays planar");
    geomlib_trimesh_free(&flat);

    geomlib_trimesh noisy=make_grid();
    for (size_t i=0;i<noisy.vertex_count;++i)
        noisy.vertices[i].z=0.08f*sinf((float)(i*13));
    float before=0.0f;
    for (size_t i=0;i<noisy.vertex_count;++i) before+=noisy.vertices[i].z*noisy.vertices[i].z;
    before=std::sqrt(before/noisy.vertex_count);
    rc=geomlib_mesh_bilateral_smooth(&noisy,nullptr);
    float after=0.0f;
    for (size_t i=0;i<noisy.vertex_count;++i) after+=noisy.vertices[i].z*noisy.vertices[i].z;
    after=std::sqrt(after/noisy.vertex_count);
    std::printf("bilateral noisy rc=%d z RMS before=%.4f after=%.4f\n",rc,before,after);
    check(rc==GEOMLIB_MESH_SUCCESS,"noisy smoothing succeeds");
    check(after<before,"height noise is reduced");

    bool finite=true;
    for (size_t i=0;i<noisy.vertex_count;++i)
        if (!std::isfinite(noisy.vertices[i].x)||!std::isfinite(noisy.vertices[i].y)||
            !std::isfinite(noisy.vertices[i].z)) finite=false;
    check(finite,"all positions finite");
    for (size_t i=0;i<noisy.vertex_count && noisy.normals;++i)
        if (std::fabs(geomlib_vec3_length(noisy.normals[i])-1.0f)>1e-3f)
            check(false,"normals are unit length");
    geomlib_trimesh_free(&noisy);

    if (failures!=0) {
        std::printf("%d bilateral test(s) failed\n",failures);
        return 1;
    }
    std::printf("All bilateral tests passed\n");
    return 0;
}
