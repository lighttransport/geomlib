/*
 * ABF tests: angle distortion on a curved mesh does not increase and stays
 * finite; pinned UVs remain fixed.
 */
#include <geomlib/geomlib_abf.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    /* Smooth non-planar 5x5 grid patch. */
    geomlib_trimesh sphere = {};
    geomlib_trimesh_init(&sphere, 25, 32);
    for (int y=0;y<5;++y)
        for (int x=0;x<5;++x)
            geomlib_trimesh_add_vertex(&sphere,
                geomlib_vec3_create((float)x,(float)y,0.7f*sinf(0.7f*x)*sinf(0.7f*y)),
                nullptr,nullptr);
    for (int y=0;y<4;++y)
        for (int x=0;x<4;++x) {
            uint32_t i=(uint32_t)(y*5+x);
            geomlib_trimesh_add_triangle(&sphere,i,i+1,i+6);
            geomlib_trimesh_add_triangle(&sphere,i,i+6,i+5);
        }
    int rc = GEOMLIB_MESH_SUCCESS;

    geomlib_uv_coords init = geomlib_uv_coords_create(sphere.vertex_count,
                                                     sphere.triangle_count);
    rc = geomlib_parameterize_lscm_free(&sphere, 0, 4, &init);
    check(rc == GEOMLIB_MESH_SUCCESS, "free LSCM initialization succeeds");
    if (rc != GEOMLIB_MESH_SUCCESS) return 1;
    float e0 = geomlib_abf_angle_error(&sphere, &init);
    geomlib_vec2 pin0=init.uvs[0], pin1=init.uvs[4];

    geomlib_abf_params p = geomlib_abf_params_default();
    p.iterations = 40;
    geomlib_uv_coords out = {};
    rc = geomlib_parameterize_abf(&sphere, 0, 4, &p, &out);
    float e1 = rc == GEOMLIB_MESH_SUCCESS ? geomlib_abf_angle_error(&sphere,&out) : INFINITY;
    std::printf("abf rc=%d angle error before=%.6f after=%.6f tris=%zu\n",
                rc, e0, e1, out.triangle_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "ABF succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(e1 < e0, "angle distortion decreases");
        check(std::fabs(out.uvs[0].x-pin0.x)<1e-5f &&
              std::fabs(out.uvs[0].y-pin0.y)<1e-5f &&
              std::fabs(out.uvs[4].x-pin1.x)<1e-5f &&
              std::fabs(out.uvs[4].y-pin1.y)<1e-5f, "pinned UVs unchanged");
        bool finite=true;
        for (size_t i=0;i<out.vertex_count;++i)
            if (!std::isfinite(out.uvs[i].x) || !std::isfinite(out.uvs[i].y))
                finite=false;
        check(finite, "all UVs are finite");
        geomlib_uv_coords_free(&out);
    }
    geomlib_uv_coords_free(&init);
    geomlib_trimesh_free(&sphere);

    if (failures != 0) {
        std::printf("%d ABF test(s) failed\n", failures);
        return 1;
    }
    std::printf("All ABF tests passed\n");
    return 0;
}
