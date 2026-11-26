/*
 * Sparse direct solver test: SPD tridiagonal systems solved exactly.
 */
#include <geomlib/geomlib_sparse_direct.h>

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
    const size_t n = 64;
    std::vector<uint32_t> rows, cols;
    std::vector<double> vals, b, x(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        rows.push_back((uint32_t)i); cols.push_back((uint32_t)i); vals.push_back(2.0);
        if (i > 0) { rows.push_back((uint32_t)i); cols.push_back((uint32_t)(i-1)); vals.push_back(-1.0); }
        if (i+1 < n) { rows.push_back((uint32_t)i); cols.push_back((uint32_t)(i+1)); vals.push_back(-1.0); }
    }
    geomlib_spmat A;
    int rc = geomlib_spmat_from_triplets(rows.data(), cols.data(), vals.data(),
                                         vals.size(), n, n, &A);
    check(rc == GEOMLIB_SPARSE_SUCCESS, "matrix assembly");
    for (size_t i = 0; i < n; ++i) b.push_back(i < n/2 ? 1.0 : 0.0);

    size_t bw = 0;
    rc = geomlib_solve_sparse_direct(&A, b.data(), x.data(), &bw);
    double residual = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double ax = 2.0*x[i];
        if (i > 0) ax -= x[i-1];
        if (i+1 < n) ax -= x[i+1];
        residual += (ax-b[i])*(ax-b[i]);
    }
    std::printf("sparse direct rc=%d bandwidth=%zu residual=%.3e x0=%.4f\n",
                rc, bw, residual, x[0]);
    check(rc == GEOMLIB_SPARSE_SUCCESS, "direct solve succeeds");
    check(bw <= 2, "RCM bandwidth is small");
    check(residual < 1e-10, "residual is tiny");
    geomlib_spmat_free(&A);
    if (failures != 0) {
        std::printf("%d sparse direct test(s) failed\n", failures);
        return 1;
    }
    std::printf("All sparse direct tests passed\n");
    return 0;
}
