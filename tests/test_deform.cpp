/*
 * ARAP deformation test: move two corner handles of a grid and check the
 * handles reach their targets while the mesh remains finite and connected.
 */
#include <geomlib/geomlib_deform.h>
#include <geomlib/geomlib_topology.h>

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
    geomlib_trimesh_init(&grid, 25, 32);
    for (int y = 0; y < 5; ++y)
        for (int x = 0; x < 5; ++x)
            geomlib_trimesh_add_vertex(&grid,
                geomlib_vec3_create((float)x*0.25f, (float)y*0.25f, 0.0f), NULL, NULL);
    for (int y = 0; y < 4; ++y)
        for (int x = 0; x < 4; ++x) {
            uint32_t a = (uint32_t)(y*5+x), b = a+1, c = a+5, d = c+1;
            geomlib_trimesh_add_triangle(&grid, a, b, d);
            geomlib_trimesh_add_triangle(&grid, a, d, c);
        }

    geomlib_deform_handle handles[2];
    handles[0].vertex = 0;
    handles[0].target = geomlib_vec3_create(0.0f, 0.0f, 0.0f);
    handles[1].vertex = 24;
    handles[1].target = geomlib_vec3_create(1.4f, 1.3f, 0.2f);

    geomlib_trimesh out = {};
    int rc = geomlib_deform_arap(&grid, handles, 2, 8, &out);
    std::printf("arap deform rc=%d vertices=%zu\n", rc, out.vertex_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "arap deformation succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(out.vertex_count == grid.vertex_count &&
              out.triangle_count == grid.triangle_count, "topology is preserved");
        geomlib_vec3 moved = out.vertices[24];
        check(std::fabs(moved.x-1.4f) < 0.05f && std::fabs(moved.y-1.3f) < 0.05f &&
              std::fabs(moved.z-0.2f) < 0.05f, "handle reaches its target");
        bool finite = true;
        for (size_t i = 0; i < out.vertex_count; ++i)
            if (!std::isfinite(out.vertices[i].x) || !std::isfinite(out.vertices[i].y) ||
                !std::isfinite(out.vertices[i].z)) finite = false;
        check(finite, "all deformed vertices are finite");
        geomlib_trimesh_free(&out);
    }
    geomlib_trimesh_free(&grid);

    if (failures != 0) {
        std::printf("%d deformation test(s) failed\n", failures);
        return 1;
    }
    std::printf("All deformation tests passed\n");
    return 0;
}
