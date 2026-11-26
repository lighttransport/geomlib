/*
 * Mesh boolean tests: two overlapping axis-aligned boxes with known volumes.
 */
#include <geomlib/geomlib_boolean.h>
#include <geomlib/geomlib_repair.h>
#include <geomlib/geomlib_topology.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static geomlib_trimesh welded_box(float cx, float cy, float cz, float size)
{
    geomlib_trimesh raw;
    geomlib_create_box(&raw, size, size, size);
    for (size_t i = 0; i < raw.vertex_count; ++i) {
        raw.vertices[i].x += cx;
        raw.vertices[i].y += cy;
        raw.vertices[i].z += cz;
    }
    geomlib_trimesh out;
    geomlib_soup_to_mesh(raw.vertices, raw.vertex_count, raw.indices,
                         raw.triangle_count * 3, 1e-5f, &out);
    geomlib_trimesh_free(&raw);
    return out;
}

int main(void)
{
    geomlib_trimesh a = welded_box(0, 0, 0, 1.0f);
    geomlib_trimesh b = welded_box(0.5f, 0, 0, 1.0f);
    double va = geomlib_trimesh_volume(&a);
    double vb = geomlib_trimesh_volume(&b);
    std::printf("input volumes: A=%.4f B=%.4f\n", va, vb);
    check(std::fabs(fabs(va) - 1.0) < 1e-3, "box A volume is one");
    check(std::fabs(fabs(vb) - 1.0) < 1e-3, "box B volume is one");

    const geomlib_boolean_op ops[] = {
        GEOMLIB_BOOLEAN_UNION, GEOMLIB_BOOLEAN_INTERSECTION, GEOMLIB_BOOLEAN_DIFFERENCE
    };
    const double expected[] = {1.5, 0.5, 0.5};
    const char *names[] = {"union", "intersection", "difference"};

    for (int k = 0; k < 3; ++k) {
        geomlib_trimesh out;
        int rc = geomlib_mesh_boolean(&a, &b, ops[k], 0.0f, &out);
        check(rc == GEOMLIB_MESH_SUCCESS, "boolean succeeds");
        if (rc != GEOMLIB_MESH_SUCCESS) continue;
        double vol = fabs(geomlib_trimesh_volume(&out));
        std::printf("%s: triangles=%zu volume=%.4f expected=%.4f\n",
                    names[k], out.triangle_count, vol, expected[k]);
        check(out.triangle_count > 0, "boolean output has triangles");
        check(std::fabs(vol - expected[k]) < 0.08, "boolean volume matches");
        geomlib_trimesh_free(&out);
    }

    geomlib_trimesh_free(&a);
    geomlib_trimesh_free(&b);
    if (failures != 0) {
        std::printf("%d boolean test(s) failed\n", failures);
        return 1;
    }
    std::printf("All boolean tests passed\n");
    return 0;
}
