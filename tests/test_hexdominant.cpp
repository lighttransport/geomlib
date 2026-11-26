/*
 * Hex-dominant meshing test: a tetrahedralized box yields hex cells and the
 * original tetrahedra are partitioned into hexes plus remaining tets.
 */
#include <geomlib/geomlib_hexdominant.h>

#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    /* Five tetrahedra around the shared edge (0,1) with six ring vertices:
     * the canonical 5-tet-to-hex merge configuration. */
    geomlib_tetramesh tets;
    geomlib_tetramesh_init(&tets, 8, 5);
    for (int i = 0; i < 8; ++i)
        geomlib_tetramesh_add_vertex(&tets,
            geomlib_vec3_create((float)i, (float)(i%3), (float)(i%2)));
    static const uint32_t t5[5][4] = {
        {0,1,2,3}, {0,1,3,4}, {0,1,4,5}, {0,1,5,6}, {0,1,6,7}
    };
    for (int i = 0; i < 5; ++i)
        geomlib_tetramesh_add_tetrahedron(&tets, t5[i][0], t5[i][1], t5[i][2], t5[i][3]);
    int rc = GEOMLIB_MESH_SUCCESS;

    geomlib_hexmesh hm = geomlib_hexmesh_create();
    rc = geomlib_hex_dominant(&tets, &hm);
    std::printf("hex dominant rc=%d hexes=%zu remaining_tets=%zu (input %zu)\n",
                rc, hm.hex_count, hm.tet_count, tets.tetrahedron_count);
    check(rc == GEOMLIB_MESH_SUCCESS, "hex dominant conversion succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        bool valid = true;
        for (size_t h = 0; h < hm.hex_count; ++h) {
            for (int a = 0; a < 8 && valid; ++a) {
                if (hm.hexes[h*8+a] >= hm.vertex_count) valid = false;
                for (int b = a+1; b < 8; ++b)
                    if (hm.hexes[h*8+a] == hm.hexes[h*8+b]) valid = false;
            }
        }
        check(valid, "hexes have eight distinct in-range vertices");
        check(hm.hex_count >= 1, "at least one hex cell found");
        check(5*hm.hex_count + hm.tet_count == tets.tetrahedron_count,
              "tetrahedra are partitioned into hexes and remaining tets");
        geomlib_hexmesh_free(&hm);
    }
    geomlib_tetramesh_free(&tets);
    if (failures != 0) {
        std::printf("%d hex dominant test(s) failed\n", failures);
        return 1;
    }
    std::printf("All hex dominant tests passed\n");
    return 0;
}
