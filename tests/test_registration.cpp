/*
 * Rigid ICP registration test: rotate and translate a sphere, recover the
 * transform and check alignment error.
 */
#include <geomlib/geomlib_registration.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh target;
    geomlib_create_sphere(&target, 1.0f, 16, 8);
    geomlib_trimesh source;
    geomlib_create_sphere(&source, 1.0f, 16, 8);

    /* Perturb source by a known rotation + translation. */
    geomlib_mat4 perturb = geomlib_mat4_identity();
    {
        float a = 0.4f;
        perturb.m[0] = cosf(a); perturb.m[1] = 0; perturb.m[2] = -sinf(a);
        perturb.m[4] = 0; perturb.m[5] = 1; perturb.m[6] = 0;
        perturb.m[8] = sinf(a); perturb.m[9] = 0; perturb.m[10] = cosf(a);
        perturb.m[12] = 0.3f; perturb.m[13] = -0.2f; perturb.m[14] = 0.1f;
    }
    geomlib_transform_mesh(&source, &perturb);

    geomlib_icp_params params = geomlib_icp_params_default();
    params.max_iterations = 40;
    geomlib_mat4 transform;
    float rms = geomlib_icp_rigid(&source, &target, &params, &transform);
    std::printf("initial-versus-final RMS: %.6f\n", rms);
    check(std::isfinite(rms), "ICP returns a finite RMS");
    check(rms < 0.05f, "ICP aligns the sphere within 0.05");

    geomlib_transform_mesh(&source, &transform);
    float final_rms = 0.0f;
    geomlib_bvh bvh;
    geomlib_bvh_build(&bvh, &target, GEOMLIB_BVH_QUALITY_MEDIUM);
    for (size_t i = 0; i < source.vertex_count; ++i) {
        geomlib_vec3 c; float d;
        if (geomlib_bvh_closest_point(&bvh, source.vertices[i], &c, &d, NULL))
            final_rms += d*d;
    }
    final_rms = sqrtf(final_rms / source.vertex_count);
    std::printf("post-transform RMS: %.6f\n", final_rms);
    check(final_rms < 0.05f, "returned transform aligns the source");
    geomlib_bvh_free(&bvh);

    geomlib_trimesh_free(&source);
    geomlib_trimesh_free(&target);
    if (failures != 0) {
        std::printf("%d registration test(s) failed\n", failures);
        return 1;
    }
    std::printf("All registration tests passed\n");
    return 0;
}
