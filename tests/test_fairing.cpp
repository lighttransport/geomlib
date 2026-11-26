/*
 * Implicit fairing test: noisy grid vertices become smoother while boundary
 * vertices stay fixed.
 */
#include <geomlib/geomlib_fairing.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh grid;
    geomlib_trimesh_init(&grid, 36, 50);
    for (int y = 0; y < 6; ++y)
        for (int x = 0; x < 6; ++x) {
            float z = ((x*7+y*13)%5 - 2) * 0.02f;
            if (x==0 || y==0 || x==5 || y==5) z = 0.0f;
            geomlib_trimesh_add_vertex(&grid,
                geomlib_vec3_create((float)x*0.2f, (float)y*0.2f, z), NULL, NULL);
        }
    for (int y = 0; y < 5; ++y)
        for (int x = 0; x < 5; ++x) {
            uint32_t a = (uint32_t)(y*6+x), b=a+1, c=a+6, d=c+1;
            geomlib_trimesh_add_triangle(&grid, a, b, d);
            geomlib_trimesh_add_triangle(&grid, a, d, c);
        }
    double before = 0.0;
    for (size_t i = 0; i < grid.vertex_count; ++i) before += std::fabs(grid.vertices[i].z);

    geomlib_fairing_params params = geomlib_fairing_params_default();
    params.lambda = 0.05f;
    params.iterations = 3;
    int rc = geomlib_fair_mesh(&grid, &params);
    std::printf("fairing rc=%d mean|z| %.5f -> ", rc, before/grid.vertex_count);
    double after = 0.0;
    for (size_t i = 0; i < grid.vertex_count; ++i) after += std::fabs(grid.vertices[i].z);
    std::printf("%.5f\n", after/grid.vertex_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "implicit fairing succeeds");
    check(after < before, "fairing reduces noise");

    bool boundary_fixed = true;
    for (int x = 0; x < 6; ++x) {
        for (int y = 0; y < 6; ++y) {
            if (x && y && x!=5 && y!=5) continue;
            size_t i = (size_t)(y*6+x);
            if (std::fabs(grid.vertices[i].z) > 1e-5f) boundary_fixed = false;
        }
    }
    check(boundary_fixed, "boundary vertices remain fixed");
    geomlib_trimesh_free(&grid);
    if (failures != 0) {
        std::printf("%d fairing test(s) failed\n", failures);
        return 1;
    }
    std::printf("All fairing tests passed\n");
    return 0;
}
