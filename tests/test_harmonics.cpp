/*
 * Manifold harmonics tests: spectrum ordering, zero eigenvalue for the
 * constant function, and the M-orthogonality of the computed eigenfunctions.
 */
#include <geomlib/geomlib_harmonics.h>
#include <geomlib/geomlib_sparse.h>

#include <cmath>
#include <cstdio>

static int failures = 0;
static void check(bool ok, const char *name)
{
    if (!ok) { std::printf("FAIL %s\n", name); failures++; }
}

int main(void)
{
    geomlib_trimesh sphere;
    geomlib_create_sphere(&sphere, 1.0f, 16, 8);
    geomlib_spmat mass;
    geomlib_massmatrix(&sphere, &mass);

    geomlib_harmonics_result r = {0};
    int rc = geomlib_compute_manifold_harmonics(&sphere, 6, &r);
    std::printf("harmonics rc=%d count=%zu lambda0=%.6f lambda1=%.6f\n",
                rc, r.count, rc == 0 ? r.eigenvalues[0] : -1.0,
                rc == 0 ? r.eigenvalues[1] : -1.0);
    check(rc == GEOMLIB_MESH_SUCCESS, "harmonics computation succeeds");
    if (rc == GEOMLIB_MESH_SUCCESS) {
        check(r.count == 6 && r.vertex_count == sphere.vertex_count, "result dimensions");
        check(std::fabs(r.eigenvalues[0]) < 1e-4, "first eigenvalue is zero");
        for (size_t i = 1; i < r.count; ++i)
            check(r.eigenvalues[i] >= r.eigenvalues[i-1] - 1e-6,
                  "eigenvalues are ascending");
        check(r.eigenvalues[5] > 0.0, "higher modes have positive eigenvalues");

        /* M-orthogonality: <phi_i, phi_j>_M ~ 0 for i != j. */
        if (r.count >= 2) {
            double dot = 0.0;
            for (size_t v = 0; v < sphere.vertex_count; ++v) {
                double m = 0.0;
                for (size_t k = mass.rowptr[v]; k < mass.rowptr[v+1]; ++k)
                    if (mass.colidx[k] == v) m = mass.values[k];
                dot += m * r.eigenvectors[v] * r.eigenvectors[sphere.vertex_count+v];
            }
            check(std::fabs(dot) < 1e-5, "distinct modes are mass-orthogonal");
        }
        geomlib_harmonics_result_free(&r);
    }
    geomlib_spmat_free(&mass);
    geomlib_trimesh_free(&sphere);
    if (failures != 0) {
        std::printf("%d harmonics test(s) failed\n", failures);
        return 1;
    }
    std::printf("All harmonics tests passed\n");
    return 0;
}
