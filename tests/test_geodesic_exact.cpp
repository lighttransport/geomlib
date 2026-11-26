/*
 * Continuous geodesic unfolding tests.
 */
#include <geomlib/geomlib_geodesic_exact.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    /* Unit square split along the (0,2) diagonal: path 0->1->2 unfolds to
     * the exact diagonal length sqrt(2). */
    geomlib_trimesh square;
    geomlib_trimesh_init(&square, 4, 2);
    geomlib_trimesh_add_vertex(&square, geomlib_vec3_create(0,0,0),0,0);
    geomlib_trimesh_add_vertex(&square, geomlib_vec3_create(1,0,0),0,0);
    geomlib_trimesh_add_vertex(&square, geomlib_vec3_create(1,1,0),0,0);
    geomlib_trimesh_add_vertex(&square, geomlib_vec3_create(0,1,0),0,0);
    geomlib_trimesh_add_triangle(&square, 0,1,2);
    geomlib_trimesh_add_triangle(&square, 0,2,3);

    double len = 0.0;
    int rc = geomlib_geodesic_length_unfolding(&square, 0, 2, &len);
    std::printf("square diagonal rc=%d length=%.5f expected=%.5f\n",
                rc, len, std::sqrt(2.0));
    check(rc == GEOMLIB_MESH_SUCCESS, "unfolding succeeds");
    check(std::fabs(len - std::sqrt(2.0)) < 1e-4, "unfolded diagonal is exact");
    geomlib_trimesh_free(&square);

    /* Sphere: geodesic between opposite-ish vertices is longer than the
     * straight chord and finite. */
    geomlib_trimesh sphere;
    geomlib_create_sphere(&sphere, 1.0f, 16, 8);
    /* Adjacent ring vertex: a two-point strip whose unfolding is the edge. */
    rc = geomlib_geodesic_length_unfolding(&sphere, 0, 17, &len);
    geomlib_vec3 a = sphere.vertices[0], b = sphere.vertices[17];
    double chord = geomlib_vec3_length(geomlib_vec3_sub(a,b));
    std::printf("sphere rc=%d length=%.4f chord=%.4f\n", rc, len, chord);
    check(rc == GEOMLIB_MESH_SUCCESS, "sphere unfolding succeeds");
    check(len >= chord - 1e-6 && std::isfinite(len), "geodesic >= chord and finite");
    geomlib_trimesh_free(&sphere);

    if (failures != 0) {
        std::printf("%d geodesic unfolding test(s) failed\n", failures);
        return 1;
    }
    std::printf("All geodesic unfolding tests passed\n");
    return 0;
}
