/*
 * Non-rigid ICP: recover a smooth deformation applied to a sphere.
 */
#include <geomlib/geomlib_nonrigid.h>
#include <geomlib/geomlib_bvh.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh src = {}, target = {}, fitted = {};
    int rc = geomlib_create_sphere(&src, 1.0f, 16, 12);
    check(rc == GEOMLIB_MESH_SUCCESS, "sphere created");
    if (rc != GEOMLIB_MESH_SUCCESS) return 1;

    geomlib_trimesh_init(&target, src.vertex_count, src.triangle_count);
    for (size_t i = 0; i < src.vertex_count; ++i) {
        geomlib_vec3 p = src.vertices[i];
        p.y += 0.30f*sinf(p.x);
        p.z += 0.15f*p.x;
        geomlib_trimesh_add_vertex(&target, p, NULL, NULL);
    }
    for (size_t i = 0; i < src.triangle_count; ++i)
        geomlib_trimesh_add_triangle(&target, src.indices[i*3],
            src.indices[i*3+1], src.indices[i*3+2]);

    geomlib_nonrigid_result fit = geomlib_nonrigid_result_create();
    rc = geomlib_nonrigid_icp(&src, &target, 40, 0.85f, &fit);
    std::printf("nonrigid rc=%d iters=%zu initial=%.4f final=%.4f\n",
                rc, fit.iterations, fit.initial_rms, fit.final_rms);
    check(rc == GEOMLIB_MESH_SUCCESS, "non-rigid ICP succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(fit.initial_rms > 0.05f, "initial mismatch is nonzero");
        check(fit.final_rms < fit.initial_rms*0.5f, "correspondence error drops >50%");
        bool finite = true;
        for (size_t i = 0; i < fit.vertex_count; ++i) {
            geomlib_vec3 d = fit.displacement[i];
            if (!std::isfinite(d.x) || !std::isfinite(d.y) || !std::isfinite(d.z))
                finite = false;
        }
        check(finite, "displacements are finite");

        rc = geomlib_nonrigid_apply(&src, &fit, &fitted);
        check(rc == GEOMLIB_MESH_SUCCESS, "deformation applies");
        if (rc == GEOMLIB_MESH_SUCCESS) {
            double s = 0.0;
            for (size_t i = 0; i < src.vertex_count; ++i)
                s += geomlib_vec3_length(geomlib_vec3_sub(
                    fitted.vertices[i], target.vertices[i]));
            double mean = s/src.vertex_count;
            std::printf("nonrigid vertex mean error=%.4f\n", mean);
            check(mean < 0.15, "fitted vertices approach target");
            geomlib_trimesh_free(&fitted);
        }
        geomlib_nonrigid_result_free(&fit);
    }

    geomlib_trimesh_free(&src);
    geomlib_trimesh_free(&target);
    if (failures != 0) {
        std::printf("%d non-rigid test(s) failed\n", failures);
        return 1;
    }
    std::printf("All non-rigid tests passed\n");
    return 0;
}
