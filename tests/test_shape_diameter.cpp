/*
 * Shape diameter test: sphere values are finite and roughly constant.
 */
#include <geomlib/geomlib_shape_diameter.h>

#include <algorithm>
#include <cfloat>
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
    geomlib_shape_diameter d = geomlib_shape_diameter_create(0);
    int rc = geomlib_compute_shape_diameter(&sphere, 12, &d);
    float mn=FLT_MAX, mx=-FLT_MAX, avg=0.0f;
    size_t used = 0;
    for (size_t i = 0; i < d.count; ++i) {
        if (!(d.values[i] > 0.0f)) continue;
        mn = std::min(mn, d.values[i]);
        mx = std::max(mx, d.values[i]);
        avg += d.values[i];
        used++;
    }
    avg /= used ? used : 1;
    std::printf("shape diameter rc=%d count=%zu min=%.4f max=%.4f avg=%.4f\n",
                rc, d.count, mn, mx, avg);
    check(rc == GEOMLIB_MESH_SUCCESS, "shape diameter succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(d.count == sphere.vertex_count, "one value per vertex");
        bool finite = true;
        for (size_t i = 0; i < d.count; ++i)
            if (!std::isfinite(d.values[i]) || d.values[i] < 0.0f) finite = false;
        check(finite, "values are finite and nonnegative");
        check(mx > 0.0f && (mx-mn) < 0.2f, "sphere values are roughly uniform");
    }
    geomlib_shape_diameter_free(&d);
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d shape diameter test(s) failed\n", failures);
        return 1;
    }
    std::printf("All shape diameter tests passed\n");
    return 0;
}
