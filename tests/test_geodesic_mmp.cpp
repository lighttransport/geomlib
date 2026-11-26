/*
 * MMP geodesic: planar zigzag strip and opposite cube corners.
 */
#include <geomlib/geomlib_geodesic_mmp.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static void vertex(geomlib_trimesh *m, float x, float y, float z)
{
    geomlib_trimesh_add_vertex(m, geomlib_vec3_create(x, y, z), nullptr, nullptr);
}

int main(void)
{
    geomlib_trimesh strip;
    double len = 0.0;

    /* Flat 2x1 strip with four triangles; the straight line from (0,0) to
     * (2,1) crosses face interiors and is shorter than any vertex path. */
    geomlib_trimesh_init(&strip, 6, 4);
    vertex(&strip, 0, 0, 0);
    vertex(&strip, 1, 0, 0);
    vertex(&strip, 2, 0, 0);
    vertex(&strip, 0, 1, 0);
    vertex(&strip, 1, 1, 0);
    vertex(&strip, 2, 1, 0);
    geomlib_trimesh_add_triangle(&strip, 0, 1, 4);
    geomlib_trimesh_add_triangle(&strip, 0, 4, 3);
    geomlib_trimesh_add_triangle(&strip, 1, 2, 5);
    geomlib_trimesh_add_triangle(&strip, 1, 5, 4);

    int rc = geomlib_geodesic_mmp(&strip, 0, 5, &len);
    std::printf("mmp strip rc=%d len=%.6f expect=%.6f\n",
                rc, len, std::sqrt(5.0));
    check(rc == GEOMLIB_MESH_SUCCESS && std::fabs(len - std::sqrt(5.0)) < 1e-4,
          "planar strip gives the Euclidean straight line");
    geomlib_trimesh_free(&strip);

    /* Unit cube, opposite corners: shortest surface path is sqrt(5). */
    geomlib_trimesh cube;
    geomlib_trimesh_init(&cube, 8, 12);
    vertex(&cube, 0, 0, 0); vertex(&cube, 1, 0, 0);
    vertex(&cube, 1, 1, 0); vertex(&cube, 0, 1, 0);
    vertex(&cube, 0, 0, 1); vertex(&cube, 1, 0, 1);
    vertex(&cube, 1, 1, 1); vertex(&cube, 0, 1, 1);
    geomlib_trimesh_add_triangle(&cube, 0, 1, 2);
    geomlib_trimesh_add_triangle(&cube, 0, 2, 3);
    geomlib_trimesh_add_triangle(&cube, 4, 6, 5);
    geomlib_trimesh_add_triangle(&cube, 4, 7, 6);
    geomlib_trimesh_add_triangle(&cube, 0, 4, 5);
    geomlib_trimesh_add_triangle(&cube, 0, 5, 1);
    geomlib_trimesh_add_triangle(&cube, 3, 2, 6);
    geomlib_trimesh_add_triangle(&cube, 3, 6, 7);
    geomlib_trimesh_add_triangle(&cube, 0, 3, 7);
    geomlib_trimesh_add_triangle(&cube, 0, 7, 4);
    geomlib_trimesh_add_triangle(&cube, 1, 5, 6);
    geomlib_trimesh_add_triangle(&cube, 1, 6, 2);

    rc = geomlib_geodesic_mmp(&cube, 0, 6, &len);
    std::printf("mmp cube rc=%d len=%.6f expect=%.6f\n",
                rc, len, std::sqrt(5.0));
    check(rc == GEOMLIB_MESH_SUCCESS && std::fabs(len - std::sqrt(5.0)) < 1e-4,
          "opposite cube corners give sqrt(5)");
    check(geomlib_geodesic_mmp(&cube, 0, 0, &len) == GEOMLIB_MESH_SUCCESS &&
          len == 0.0, "zero distance for identical vertices");
    geomlib_trimesh_free(&cube);

    if (failures != 0) {
        std::printf("%d MMP geodesic test(s) failed\n", failures);
        return 1;
    }
    std::printf("All MMP geodesic tests passed\n");
    return 0;
}
