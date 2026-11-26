/*
 * Robust topology repair tests: polygon-soup reconstruction and non-manifold
 * edge detachment.
 */
#include <geomlib/geomlib_repair.h>
#include <geomlib/geomlib_topology.h>

#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    /* Soup of a tetrahedron with duplicated points, reversed faces and a
     * duplicate triangle. */
    {
        geomlib_vec3 p[8] = {
            {0,0,0},{1,0,0},{0,1,0},{0,0,1},
            {0,0,0},{1,0,0},{0,1,0},{0,0,1}
        };
        uint32_t idx[] = {
            0,1,2, 0,2,3, 0,3,1, 1,3,2,
            4,6,5, 0,1,2
        };
        geomlib_trimesh mesh;
        int rc = geomlib_soup_to_mesh(p, 8, idx, 18, 1e-5f, &mesh);
        check(rc == GEOMLIB_MESH_SUCCESS, "soup to mesh rc");
        if (rc == GEOMLIB_MESH_SUCCESS) {
            check(mesh.vertex_count == 4, "duplicate points welded");
            check(mesh.triangle_count == 4, "duplicate face removed");
            geomlib_manifold_analysis a = geomlib_analyze_manifold(&mesh);
            check(a.is_closed && a.is_manifold && a.is_oriented,
                  "reconstructed soup is closed manifold with consistent winding");
            geomlib_manifold_analysis_free(&a);
            geomlib_trimesh_free(&mesh);
        }
    }

    /* Two tetrahedra sharing an edge: that edge has four incident faces. */
    {
        geomlib_trimesh mesh;
        geomlib_trimesh_init(&mesh, 6, 8);
        geomlib_vec3 v[6] = {
            {0,0,0},{1,0,0},{0,1,0},{0,0,1},{0,-1,0},{0,0,-1}
        };
        for (int i = 0; i < 6; ++i)
            geomlib_trimesh_add_vertex(&mesh, v[i], 0, 0);
        /* Shared edge is (0,1). */
        geomlib_trimesh_add_triangle(&mesh, 0,1,2);
        geomlib_trimesh_add_triangle(&mesh, 0,1,3);
        geomlib_trimesh_add_triangle(&mesh, 0,1,4);
        geomlib_trimesh_add_triangle(&mesh, 0,1,5);
        geomlib_trimesh_add_triangle(&mesh, 2,3,4);
        geomlib_trimesh_add_triangle(&mesh, 2,4,5);
        geomlib_trimesh_add_triangle(&mesh, 3,5,4);
        geomlib_trimesh_add_triangle(&mesh, 3,2,5);

        size_t detached = 0;
        check(geomlib_repair_nonmanifold_edges(&mesh, &detached) == GEOMLIB_MESH_SUCCESS,
              "edge repair rc");
        check(detached >= 2, "over-used edge was detached");
        /* After repair no edge should have more than two incident faces. */
        geomlib_manifold_analysis a = geomlib_analyze_manifold(&mesh);
        check(a.nonmanifold_count == 0, "no non-manifold elements remain");
        geomlib_manifold_analysis_free(&a);
        geomlib_trimesh_free(&mesh);
    }

    if (failures != 0) {
        std::printf("%d repair test(s) failed\n", failures);
        return 1;
    }
    std::printf("All repair tests passed\n");
    return 0;
}
