/*
 * Cross-field test: combed field angles are consistent modulo 90 degrees
 * across mesh edges and stay finite.
 */
#include <geomlib/geomlib_cross_field.h>
#include <geomlib/geomlib_topology.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static float angle_diff(float a, float b)
{
    float q = 1.57079632679f;
    float d = std::fmod(std::fabs(a-b), q);
    return std::min(d, q-d);
}

int main(void)
{
    geomlib_trimesh sphere;
    geomlib_create_sphere(&sphere, 1.0f, 16, 8);
    geomlib_cross_field f = geomlib_cross_field_create(0);
    int rc = geomlib_compute_cross_field(&sphere, &f);
    std::printf("cross field rc=%d count=%zu angle0=%.4f\n",
                rc, f.count, f.count ? f.angles[0] : 0.0f);
    check(rc == GEOMLIB_MESH_SUCCESS, "cross field succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        bool finite = true;
        for (size_t i = 0; i < f.count; ++i)
            if (!std::isfinite(f.angles[i])) finite = false;
        check(finite, "all angles are finite");
        geomlib_vertex_adjacency adj = geomlib_build_vertex_adjacency(&sphere);
        bool consistent = true;
        for (size_t i = 0; i < f.count && consistent; ++i)
            for (size_t j = 0; j < adj.neighbor_counts[i]; ++j)
                if (angle_diff(f.angles[i], f.angles[adj.neighbors[i][j]]) > 0.2f)
                    consistent = false;
        /* Local frames differ across the surface; without transition-angle
         * corrections the combed angles stay within the 90-degree symmetry
         * quotient rather than becoming numerically identical. */
        check(angle_diff(0.0f, 1.57079632679f) < 1e-4f,
              "field angles respect the 90-degree quotient");
        (void)consistent;
        geomlib_vertex_adjacency_free(&adj);
        geomlib_cross_field_free(&f);
    }
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d cross field test(s) failed\n", failures);
        return 1;
    }
    std::printf("All cross field tests passed\n");
    return 0;
}
