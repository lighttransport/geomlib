/*
 * Voxelization tests: a unit box gives the expected volume and a thin
 * surface shell; a sphere's volume approximates 4*pi*r^3/3.
 */
#include <geomlib/geomlib_voxelize.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh box = {};
    int rc = geomlib_create_box(&box, 1.0f, 1.0f, 1.0f);
    check(rc == GEOMLIB_MESH_SUCCESS, "box created");
    if (rc != GEOMLIB_MESH_SUCCESS) return 1;

    geomlib_voxel_grid g = geomlib_voxelize_mesh(&box, 24, 24, 24);
    check(g.surface != nullptr && g.volume != nullptr, "voxel grid created");
    double vol = (double)g.volume_count*g.voxel_size*g.voxel_size*g.voxel_size;
    std::printf("box voxelize surface=%zu volume=%zu voxel=%.5f volume_est=%.4f\n",
                g.surface_count, g.volume_count, g.voxel_size, vol);
    check(g.surface_count > 100, "surface shell is populated");
    check(std::fabs(vol-1.0) < 0.12, "box volume is approximately correct");
    /* Surface voxels should be much fewer than interior voxels for a box. */
    check(g.surface_count < g.volume_count, "surface shell thinner than volume");

    /* Every volume voxel centre must lie inside the box bounds. */
    int outside = 0;
    for (uint32_t z = 0; z < g.nz; ++z)
        for (uint32_t y = 0; y < g.ny; ++y)
            for (uint32_t x = 0; x < g.nx; ++x) {
                size_t id = geomlib_voxel_grid_index(&g, x, y, z);
                if (!g.volume[id]) continue;
                geomlib_vec3 p = geomlib_vec3_create(
                    g.origin.x+((float)x+0.5f)*g.voxel_size,
                    g.origin.y+((float)y+0.5f)*g.voxel_size,
                    g.origin.z+((float)z+0.5f)*g.voxel_size);
                if (p.x < -0.55f || p.x > 0.55f || p.y < -0.55f || p.y > 0.55f ||
                    p.z < -0.55f || p.z > 0.55f) outside++;
            }
    check(outside == 0, "interior voxels lie within the box");

    geomlib_voxel_grid_free(&g);

    geomlib_trimesh sphere = {};
    rc = geomlib_create_sphere(&sphere, 0.5f, 16, 12);
    if (rc == GEOMLIB_MESH_SUCCESS) {
        g = geomlib_voxelize_mesh(&sphere, 24, 24, 24);
        vol = (double)g.volume_count*g.voxel_size*g.voxel_size*g.voxel_size;
        double expect = 4.0*M_PI*0.125/3.0;
        std::printf("sphere voxelize surface=%zu volume=%zu volume_est=%.4f expect=%.4f\n",
                    g.surface_count, g.volume_count, vol, expect);
        check(std::fabs(vol-expect) < 0.12, "sphere volume is approximately correct");
        geomlib_voxel_grid_free(&g);
        geomlib_trimesh_free(&sphere);
    }
    geomlib_trimesh_free(&box);
    if (failures != 0) {
        std::printf("%d voxelization test(s) failed\n", failures);
        return 1;
    }
    std::printf("All voxelization tests passed\n");
    return 0;
}
