/*
 * Broad-phase overlap pairs: identical boxes overlap, distant boxes do not,
 * and every reported pair has valid triangle indices.
 */
#include <geomlib/geomlib_mesh_overlap.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static void shift(geomlib_trimesh *m, float x, float y, float z)
{
    for (size_t i = 0; i < m->vertex_count; ++i) {
        m->vertices[i].x += x;
        m->vertices[i].y += y;
        m->vertices[i].z += z;
    }
}

int main(void)
{
    geomlib_trimesh a, b;
    geomlib_create_box(&a, 1.0f, 1.0f, 1.0f);
    geomlib_create_box(&b, 1.0f, 1.0f, 1.0f);

    /* Identically placed boxes: many triangles have overlapping AABBs. */
    geomlib_overlap_pairs p = geomlib_overlap_pairs_create();
    int rc = geomlib_aabb_overlap_pairs(&a, &a, &p);
    std::printf("same-box rc=%d pairs=%zu (12 tris each)\n", rc, p.count);
    check(rc == GEOMLIB_MESH_SUCCESS && p.count > 0, "identical boxes overlap");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        bool valid = true;
        for (size_t i = 0; i < p.count; ++i) {
            if (p.pairs[i * 2] >= a.triangle_count ||
                p.pairs[i * 2 + 1] >= a.triangle_count) valid = false;
        }
        check(valid, "reported indices are in range");
    }
    geomlib_overlap_pairs_free(&p);

    /* Separated box: no overlap. */
    shift(&b, 5.0f, 0.0f, 0.0f);
    rc = geomlib_aabb_overlap_pairs(&a, &b, &p);
    std::printf("far-box rc=%d pairs=%zu\n", rc, p.count);
    check(rc == GEOMLIB_MESH_SUCCESS && p.count == 0,
          "distant boxes do not overlap");
    geomlib_overlap_pairs_free(&p);

    /* Overlapping boxes: some overlap, all indices in range. */
    geomlib_trimesh_free(&b);
    geomlib_create_box(&b, 1.0f, 1.0f, 1.0f);
    shift(&b, 0.5f, 0.0f, 0.0f);
    rc = geomlib_aabb_overlap_pairs(&a, &b, &p);
    std::printf("near-box rc=%d pairs=%zu\n", rc, p.count);
    check(rc == GEOMLIB_MESH_SUCCESS && p.count > 0,
          "overlapping boxes report pairs");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        bool valid = true;
        for (size_t i = 0; i < p.count; ++i) {
            if (p.pairs[i * 2] >= a.triangle_count ||
                p.pairs[i * 2 + 1] >= b.triangle_count) valid = false;
        }
        check(valid, "reported cross-mesh indices are in range");
    }
    geomlib_overlap_pairs_free(&p);

    geomlib_trimesh_free(&a);
    geomlib_trimesh_free(&b);
    if (failures != 0) {
        std::printf("%d mesh overlap test(s) failed\n", failures);
        return 1;
    }
    std::printf("All mesh overlap tests passed\n");
    return 0;
}
