/*
 * Octree tests: radius and kNN queries agree with brute force.
 */
#include <geomlib/geomlib_octree.h>

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
    for (int i = 0; i < 240; ++i) {
        float u = i*0.618033f*6.283185f;
        float r = (i%19)/18.0f;
        pts.push_back(geomlib_vec3_create(r*cosf(u), r*sinf(u), ((i*5)%13)/12.0f));
    }
    geomlib_octree *t = geomlib_octree_build(pts.data(), pts.size(), NULL);
    check(t != NULL, "octree builds");
    if (!t) return 1;
    geomlib_vec3 q = geomlib_vec3_create(0.2f, 0.1f, 0.4f);

    uint32_t idx[6];
    float dst[6];
    size_t n = geomlib_octree_knn(t, q, 6, idx, dst);
    float best = INFINITY;
    for (auto &p : pts)
        best = std::min(best, geomlib_vec3_length(geomlib_vec3_sub(p,q)));
    std::printf("knn count=%zu nearest=%.4f brute=%.4f\n", n, dst[0], best);
    check(n == 6, "six neighbours returned");
    check(std::fabs(dst[0]-best) < 1e-4f, "nearest neighbour matches brute force");
    bool sorted = true;
    for (size_t i = 1; i < n; ++i) if (dst[i] < dst[i-1]) sorted = false;
    check(sorted, "neighbours sorted by distance");

    uint32_t ridx[512];
    float rdst[512];
    size_t rn = geomlib_octree_radius(t, q, 0.6f, ridx, rdst, 512);
    size_t brute = 0;
    for (auto &p : pts)
        if (geomlib_vec3_length(geomlib_vec3_sub(p,q)) <= 0.6f) brute++;
    std::printf("radius count=%zu brute=%zu\n", rn, brute);
    check(rn == brute, "radius count matches brute force");

    geomlib_octree_free(t);
    if (failures != 0) {
        std::printf("%d octree test(s) failed\n", failures);
        return 1;
    }
    std::printf("All octree tests passed\n");
    return 0;
}
