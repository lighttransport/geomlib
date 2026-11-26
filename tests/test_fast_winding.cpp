/*
 * Fast winding number tests: hierarchical approximation should agree with
 * the exact generalized winding number on a sphere.
 */
#include <geomlib/geomlib_fast_winding.h>
#include <geomlib/geomlib_sdf.h>

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
    geomlib_winding_tree *tree = geomlib_winding_tree_build(&sphere, 24);
    check(tree != NULL, "winding tree builds");
    if (tree) {
        geomlib_vec3 inside = geomlib_vec3_create(0.2f, 0.1f, 0.1f);
        geomlib_vec3 outside = geomlib_vec3_create(2.0f, 1.0f, 0.5f);
        float fi = geomlib_fast_winding_number(tree, inside);
        float fo = geomlib_fast_winding_number(tree, outside);
        float ei = geomlib_compute_winding_number(&sphere, inside);
        float eo = geomlib_compute_winding_number(&sphere, outside);
        std::printf("fast inside=%.4f exact=%.4f outside=%.4f exact=%.4f\n",
                    fi, ei, fo, eo);
        check(std::fabs(fi-ei) < 0.05f, "inside winding matches exact");
        check(std::fabs(fo-eo) < 0.05f, "outside winding matches exact");
        check(std::fabs(fi) > 0.9f && std::fabs(fo) < 0.05f,
              "inside/outside classification");
        geomlib_winding_tree_free(tree);
    }
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d fast winding test(s) failed\n", failures);
        return 1;
    }
    std::printf("All fast winding tests passed\n");
    return 0;
}
