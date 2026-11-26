/*
 * Blue-noise surface sampling test: separation, surface proximity and
 * determinism with a fixed RNG seed.
 */
#include <geomlib/geomlib_bluenoise.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh sphere;
    geomlib_create_sphere(&sphere, 1.0f, 24, 12);
    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 7);
    geomlib_sample_set samples = {0};
    float min_distance = 0.12f;
    size_t n = geomlib_blue_noise_surface(&sphere, min_distance, &rng, &samples);
    std::printf("samples=%zu min_distance=%.3f\n", n, min_distance);
    check(n > 50, "generated a useful number of samples");

    int min_ok = 1, on_surface = 1;
    for (size_t i = 0; i < samples.count; ++i) {
        geomlib_vec3 p = samples.points[i];
        check(std::fabs(std::sqrt(p.x*p.x+p.y*p.y+p.z*p.z) - 1.0f) < 0.05f,
              "samples lie near the unit sphere");
        for (size_t j = i+1; j < samples.count; ++j) {
            geomlib_vec3 d = geomlib_vec3_sub(p, samples.points[j]);
            if (geomlib_vec3_length(d) < min_distance * 0.999f) min_ok = 0;
        }
    }
    (void)on_surface;
    check(min_ok, "Poisson-disk separation is respected");

    geomlib_pcg32 rng2;
    geomlib_pcg32_init(&rng2, 42, 7);
    geomlib_sample_set again = {0};
    size_t n2 = geomlib_blue_noise_surface(&sphere, min_distance, &rng2, &again);
    check(n2 == n, "sampling is deterministic for a fixed seed");
    if (n == n2 && n > 0) {
        bool same = true;
        for (size_t i = 0; i < n; ++i) {
            geomlib_vec3 a = samples.points[i], b = again.points[i];
            if (fabsf(a.x-b.x) > 1e-6f || fabsf(a.y-b.y) > 1e-6f || fabsf(a.z-b.z) > 1e-6f)
                same = false;
        }
        check(same, "sampled points match exactly for the same seed");
    }

    geomlib_sample_set_free(&samples);
    geomlib_sample_set_free(&again);
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d blue-noise test(s) failed\n", failures);
        return 1;
    }
    std::printf("All blue-noise tests passed\n");
    return 0;
}
