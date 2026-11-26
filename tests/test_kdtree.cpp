/*
 * k-d tree tests: kNN and radius queries agree with brute force.
 */
#include <geomlib/geomlib_kdtree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    std::vector<geomlib_vec3> pts;
    for (int i = 0; i < 200; ++i) {
        float u = i*0.618033f*6.283185f;
        float r = (i%17)/16.0f;
        pts.push_back(geomlib_vec3_create(r*cosf(u), r*sinf(u), ((i*7)%11)/10.0f));
    }
    geomlib_kdtree *tree = geomlib_kdtree_build(pts.data(), pts.size());
    check(tree != NULL, "k-d tree builds");
    if (!tree) return 1;

    geomlib_vec3 q = geomlib_vec3_create(0.1f, 0.2f, 0.3f);
    uint32_t idx[5];
    float dist[5];
    size_t n = geomlib_kdtree_knn(tree, q, 5, idx, dist);
    std::printf("knn count=%zu d0=%.4f d4=%.4f\n", n, dist[0], dist[4]);
    check(n == 5, "five neighbours returned");
    bool sorted = true;
    for (size_t i = 1; i < n; ++i) if (dist[i] < dist[i-1]) sorted = false;
    check(sorted, "neighbours are sorted by distance");
    /* Brute-force nearest. */
    float best = INFINITY;
    for (auto &p : pts) {
        geomlib_vec3 d = geomlib_vec3_sub(p,q);
        best = std::min(best, geomlib_vec3_length(d));
    }
    check(std::fabs(dist[0]-best) < 1e-4f, "nearest neighbour matches brute force");

    uint32_t ridx[256];
    float rdist[256];
    size_t rn = geomlib_kdtree_radius(tree, q, 0.5f, ridx, rdist, 256);
    size_t brute = 0;
    for (auto &p : pts)
        if (geomlib_vec3_length(geomlib_vec3_sub(p,q)) <= 0.5f) brute++;
    std::printf("radius count=%zu brute=%zu\n", rn, brute);
    check(rn == brute, "radius query matches brute force");

    geomlib_kdtree_free(tree);
    if (failures != 0) {
        std::printf("%d k-d tree test(s) failed\n", failures);
        return 1;
    }
    std::printf("All k-d tree tests passed\n");
    return 0;
}
