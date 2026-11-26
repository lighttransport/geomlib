/*
 * Quad remeshing tests: a regular triangulated grid should pair into quads;
 * each quad must have four distinct vertices and the element count must match.
 */
#include <geomlib/geomlib_quadremesh.h>

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
        for (int x = 0; x < 6; ++x)
            geomlib_trimesh_add_vertex(&grid,
                geomlib_vec3_create((float)x*0.2f, (float)y*0.2f, 0.0f), NULL, NULL);
    for (int y = 0; y < 5; ++y)
        for (int x = 0; x < 5; ++x) {
            uint32_t a = (uint32_t)(y*6+x), b = a+1, c = a+6, d = c+1;
            geomlib_trimesh_add_triangle(&grid, a, b, d);
            geomlib_trimesh_add_triangle(&grid, a, d, c);
        }

    geomlib_quadremesh_params params = geomlib_quadremesh_params_default();
    params.max_dihedral = 1.0f;    /* planar grid */
    params.smoothing_iterations = 2;
    geomlib_quadmesh out = {};
    int rc = geomlib_remesh_to_quads(&grid, &params, &out);
    std::printf("quadremesh rc=%d quads=%zu triangles=%zu\n",
                rc, out.quad_count, out.triangle_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "quad remeshing succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(out.quad_count >= 20, "many triangles paired into quads");
        check(2*out.quad_count + out.triangle_count == grid.triangle_count,
              "element count is conserved");
        bool valid = true;
        for (size_t q = 0; q < out.quad_count; ++q) {
            for (int a = 0; a < 4; ++a) {
                if (out.quads[q*4+a] >= out.vertex_count) valid = false;
                for (int b = a+1; b < 4; ++b)
                    if (out.quads[q*4+a] == out.quads[q*4+b]) valid = false;
            }
        }
        check(valid, "quads have four distinct in-range vertices");
        geomlib_quadmesh_free(&out);
    }
    geomlib_trimesh_free(&grid);
    if (failures != 0) {
        std::printf("%d quad remeshing test(s) failed\n", failures);
        return 1;
    }
    std::printf("All quad remeshing tests passed\n");
    return 0;
}
