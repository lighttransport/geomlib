/*
 * Spatial sorting tests: Morton/Hilbert encode-decode round trips, bijectivity,
 * and locality of the resulting point order.
 */
#include <geomlib/geomlib_spatial_sort.h>

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
    int bad = 0;
    for (uint32_t x = 0; x < 64; ++x)
        for (uint32_t y = 0; y < 64; ++y) {
            uint32_t dx, dy;
            geomlib_morton2_decode(geomlib_morton2_encode(x*1000+7, y*2003+11), &dx, &dy);
            if (dx != x*1000+7 || dy != y*2003+11) bad++;
        }
    std::printf("morton2 round-trip failures=%d\n", bad);
    check(bad == 0, "Morton 2D round trips");

    bad = 0;
    for (uint32_t x = 0; x < 32; ++x)
        for (uint32_t y = 0; y < 32; ++y)
            for (uint32_t z = 0; z < 32; ++z) {
                uint32_t dx, dy, dz;
                uint32_t a=(x<<6)+x, b=(y<<6)+y, c=(z<<6)+z;
                geomlib_morton3_decode(geomlib_morton3_encode(a,b,c), &dx, &dy, &dz);
                if (dx != a || dy != b || dz != c) bad++;
            }
    std::printf("morton3 round-trip failures=%d\n", bad);
    check(bad == 0, "Morton 3D round trips");

    std::vector<uint64_t> codes;
    bad = 0;
    for (uint32_t y = 0; y < 64; ++y)
        for (uint32_t x = 0; x < 64; ++x) {
            uint32_t dx, dy;
            uint64_t d = geomlib_hilbert2_encode(64, x, y);
            geomlib_hilbert2_decode(64, d, &dx, &dy);
            if (dx != x || dy != y) bad++;
            codes.push_back(d);
        }
    std::sort(codes.begin(), codes.end());
    size_t unique = codes.size();
    codes.erase(std::unique(codes.begin(), codes.end()), codes.end());
    std::printf("hilbert round-trip failures=%d unique=%zu/%zu\n",
                bad, codes.size(), unique);
    check(bad == 0 && codes.size() == unique, "Hilbert is a bijection on 64x64");

    std::vector<geomlib_vec3> pts;
    for (int i = 0; i < 400; ++i)
        pts.push_back(geomlib_vec3_create(
            ((i*37)%101)/100.0f, ((i*73)%97)/100.0f, ((i*29)%89)/100.0f));
    std::vector<uint32_t> order(pts.size());
    int rc = geomlib_spatial_sort3(pts.data(), pts.size(), order.data());
    check(rc == GEOMLIB_MESH_SUCCESS, "3D spatial sort succeeds");
    std::vector<uint8_t> seen(pts.size(), 0);
    bool perm = true;
    for (size_t i = 0; i < order.size(); ++i) {
        if (order[i] >= pts.size() || seen[order[i]]) perm = false;
        else seen[order[i]] = 1;
    }
    check(perm, "sort returns a permutation");

    double consecutive = 0.0, random = 0.0;
    for (size_t i = 1; i < order.size(); ++i)
        consecutive += geomlib_vec3_length(
            geomlib_vec3_sub(pts[order[i]], pts[order[i-1]]));
    consecutive /= (double)(order.size()-1);
    for (size_t i = 0; i < 400; ++i) {
        size_t a=(size_t)((i*17)%400), b=(size_t)((i*53+1)%400);
        random += geomlib_vec3_length(geomlib_vec3_sub(pts[a], pts[b]));
    }
    random /= 400.0;
    std::printf("spatial sort consecutive=%.4f random=%.4f\n", consecutive, random);
    check(consecutive < random*0.5, "sorted neighbours are spatially local");

    if (failures != 0) {
        std::printf("%d spatial sort test(s) failed\n", failures);
        return 1;
    }
    std::printf("All spatial sort tests passed\n");
    return 0;
}
