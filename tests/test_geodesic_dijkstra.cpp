/*
 * Dijkstra geodesic tests on a connected sphere: symmetry, geodesic >=
 * Euclidean, finite distances and multi-source dominance.
 */
#include <geomlib/geomlib_geodesic_dijkstra.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh sphere = {};
    int rc = geomlib_create_sphere(&sphere, 1.0f, 16, 12);
    check(rc == GEOMLIB_MESH_SUCCESS, "sphere created");
    if (rc != GEOMLIB_MESH_SUCCESS) return 1;

    geomlib_geodesic_distances d = geomlib_geodesic_create(sphere.vertex_count);
    rc = geomlib_geodesic_dijkstra_from_vertex(&sphere, 0, &d);
    float mx=0.0f;
    for (size_t i=0;i<sphere.vertex_count;++i) if (d.distances[i]>mx) mx=d.distances[i];
    std::printf("dijkstra rc=%d source=%.4f max=%.4f\n", rc, d.distances[0], mx);
    check(rc == GEOMLIB_MESH_SUCCESS, "Dijkstra succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(std::fabs(d.distances[0]) < 1e-6f, "source distance is zero");
        bool finite=true, lower=true;
        for (size_t i=0;i<sphere.vertex_count;++i) {
            if (!std::isfinite(d.distances[i]) || d.distances[i] > 100.0f)
                finite=false;
            float eu = geomlib_vec3_length(geomlib_vec3_sub(
                sphere.vertices[i], sphere.vertices[0]));
            if (d.distances[i] + 1e-4f < eu) lower=false;
        }
        check(finite, "all distances finite");
        check(lower, "graph distance >= Euclidean distance");

        uint32_t tgt=(uint32_t)(sphere.vertex_count/2);
        float d01=d.distances[tgt];
        geomlib_geodesic_distances back=geomlib_geodesic_create(sphere.vertex_count);
        int rc2=geomlib_geodesic_dijkstra_from_vertex(&sphere,tgt,&back);
        check(rc2==GEOMLIB_MESH_SUCCESS &&
              std::fabs(back.distances[0]-d01) < 1e-3f, "distance is symmetric");
        geomlib_geodesic_free(&back);

        uint32_t sources[2]={0,tgt};
        geomlib_geodesic_distances multi=geomlib_geodesic_create(sphere.vertex_count);
        rc2=geomlib_geodesic_dijkstra(&sphere,sources,2,&multi);
        bool no_worse=true;
        for (size_t i=0;i<sphere.vertex_count;++i)
            if (multi.distances[i] > d.distances[i]+1e-4f) no_worse=false;
        check(rc2==GEOMLIB_MESH_SUCCESS && no_worse,
              "multi-source distances are no larger");
        geomlib_geodesic_free(&multi);
    }
    geomlib_geodesic_free(&d);
    geomlib_trimesh_free(&sphere);

    if (failures != 0) {
        std::printf("%d Dijkstra test(s) failed\n", failures);
        return 1;
    }
    std::printf("All Dijkstra tests passed\n");
    return 0;
}
