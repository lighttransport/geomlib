/*
 * Mesh I/O tests: round-trip OBJ/OFF/PLY/STL through memory buffers,
 * including a polygon OBJ face and negative indices.
 */
#include <geomlib/geomlib_io.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

static int failures = 0;

static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

static void compare_meshes(const geomlib_trimesh *a, const geomlib_trimesh *b,
                           const char *name)
{
    check(a->vertex_count == b->vertex_count, "vertex count");
    check(a->triangle_count == b->triangle_count, "triangle count");
    if (a->vertex_count != b->vertex_count || a->triangle_count != b->triangle_count) return;
    for (size_t i = 0; i < a->vertex_count; ++i) {
        geomlib_vec3 x = a->vertices[i], y = b->vertices[i];
        check(std::fabs(x.x-y.x) < 1e-4 && std::fabs(x.y-y.y) < 1e-4 &&
              std::fabs(x.z-y.z) < 1e-4, name);
    }
    for (size_t i = 0; i < a->triangle_count * 3; ++i) {
        check(a->indices[i] == b->indices[i], name);
    }
}

int main(void)
{
    /* Create a tetrahedron. */
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 4, 4);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0,0,0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1,0,0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0,1,0), NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0,0,1), NULL, NULL);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 3);
    geomlib_trimesh_add_triangle(&mesh, 0, 2, 3);
    geomlib_trimesh_add_triangle(&mesh, 1, 2, 3);

    const char *formats[] = {"obj", "off", "ply", "plyb", "stl"};
    for (int f = 0; f < 5; ++f) {
        void *buf = NULL;
        size_t size = 0;
        int rc = geomlib_mesh_save_mem(&mesh, formats[f], &buf, &size);
        check(rc == GEOMLIB_IO_SUCCESS && buf && size > 0, "save rc");
        if (rc != GEOMLIB_IO_SUCCESS || !buf) continue;
        check(std::strcmp(geomlib_mesh_detect_format(buf, size), formats[f]) == 0,
              "detect format");
        geomlib_trimesh loaded;
        rc = geomlib_mesh_load_mem(&loaded, NULL, buf, size);
        check(rc == GEOMLIB_IO_SUCCESS, "load rc");
        if (rc == GEOMLIB_IO_SUCCESS) {
            if (std::strcmp(formats[f], "stl") == 0) {
                /* STL is unindexed: one vertex per triangle corner. */
                check(loaded.vertex_count == mesh.triangle_count * 3, "stl vertex count");
                check(loaded.triangle_count == mesh.triangle_count, "stl triangle count");
            } else {
                compare_meshes(&mesh, &loaded, "roundtrip");
            }
            geomlib_trimesh_free(&loaded);
        }
        free(buf);
    }

    /* OBJ with v/vt/vn indices, a quad face and negative indices. */
    {
        const char obj[] =
            "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv 0 0 1\n"
            "vn 0 0 1\nvt 0 0\n"
            "f 1/1/1 2/1/1 3/1/1 4/1/1\n"
            "f -4 -2 -5\n";
        geomlib_trimesh m2;
        int rc = geomlib_mesh_load_mem(&m2, "obj", obj, sizeof(obj) - 1);
        check(rc == GEOMLIB_IO_SUCCESS, "obj polygon rc");
        if (rc == GEOMLIB_IO_SUCCESS) {
            check(m2.vertex_count == 5, "obj vertices");
            check(m2.triangle_count == 3, "obj fan triangles");
            check(m2.indices[0] == 0 && m2.indices[1] == 1 && m2.indices[2] == 2,
                  "obj face 0");
            check(m2.indices[3] == 0 && m2.indices[4] == 2 && m2.indices[5] == 3,
                  "obj face 1");
            check(m2.indices[6] == 1 && m2.indices[7] == 3 && m2.indices[8] == 0,
                  "obj negative face");
            geomlib_trimesh_free(&m2);
        }
    }

    /* ASCII STL. */
    {
        const char stl[] =
            "solid test\nfacet normal 0 0 1\nouter loop\n"
            "vertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\n"
            "endloop\nendfacet\nendsolid test\n";
        geomlib_trimesh m3;
        int rc = geomlib_mesh_load_mem(&m3, "stl", stl, sizeof(stl) - 1);
        check(rc == GEOMLIB_IO_SUCCESS, "ascii stl rc");
        if (rc == GEOMLIB_IO_SUCCESS) {
            check(m3.vertex_count == 3 && m3.triangle_count == 1, "ascii stl counts");
            geomlib_trimesh_free(&m3);
        }
    }

    geomlib_trimesh_free(&mesh);
    if (failures != 0) {
        std::printf("%d I/O test(s) failed\n", failures);
        return 1;
    }
    std::printf("All mesh I/O tests passed\n");
    return 0;
}
