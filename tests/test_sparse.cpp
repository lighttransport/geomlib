/*
 * Sparse linear algebra tests: CSR assembly, CG, CGLS, dense QR, symmetric
 * eigensolver and the geometry operators.
 */
#include <geomlib/geomlib_sparse.h>

#include <cmath>
#include <cstdio>

static int failures = 0;

static void check(bool ok, const char *name)
{
    if (!ok) {
        std::printf("FAIL %s\n", name);
        failures++;
    }
}

static void check_near(double got, double expected, double tol, const char *name)
{
    if (std::fabs(got - expected) > tol) {
        std::printf("FAIL %s: got %.10g expected %.10g\n", name, got, expected);
        failures++;
    }
}

int main(void)
{
    /* Tridiagonal SPD matrix (4x4) with a duplicate entry to test coalescing. */
    {
        const uint32_t rows[] = {0,0, 1,1,1, 2,2,2, 3,3, 0};
        const uint32_t cols[] = {0,1, 0,1,2, 1,2,3, 2,3, 0};
        const double vals[]   = {2, -1, -1, 2, -1, -1, 2, -1, -1, 2, 0.5};
        geomlib_spmat a;
        int rc = geomlib_spmat_from_triplets(rows, cols, vals, 11, 4, 4, &a);
        check(rc == GEOMLIB_SPARSE_SUCCESS, "triplets rc");
        check(a.nnz == 10, "duplicate coalesced");
        check(a.rowptr[0] == 0 && a.rowptr[4] == 10, "rowptr bounds");

        double y[4];
        const double xin[4] = {1, 2, 3, 4};
        check(geomlib_spmat_mult_vec(&a, xin, y) == GEOMLIB_SPARSE_SUCCESS, "matvec rc");
        check_near(y[0], 2.5*1 - 2, 1e-12, "matvec row0");
        check_near(y[3], -3 + 2*4, 1e-12, "matvec row3");

        double b[4] = {1, 2, 3, 4};
        double x[4] = {0, 0, 0, 0};
        check(geomlib_solve_cg_sparse(&a, b, x, NULL) == GEOMLIB_SPARSE_SUCCESS, "CG rc");
        double ax[4];
        geomlib_spmat_mult_vec(&a, x, ax);
        for (int i = 0; i < 4; ++i) check_near(ax[i], b[i], 1e-8, "CG residual");

        /* Transpose. */
        geomlib_spmat at;
        check(geomlib_spmat_transpose(&a, &at) == GEOMLIB_SPARSE_SUCCESS, "transpose rc");
        check(at.nnz == a.nnz, "transpose nnz");
        geomlib_spmat_free(&at);
        geomlib_spmat_free(&a);
    }

    /* Overdetermined least squares: A = [1 0; 0 1; 1 1], b = [1, 2, 4].
     * Exact solution is [1.5, 2.5] (residual [0.5, 0.5, -0.5]... verify). */
    {
        const uint32_t rows[] = {0, 1, 2, 2};
        const uint32_t cols[] = {0, 1, 0, 1};
        const double vals[]   = {1, 1, 1, 1};
        geomlib_spmat a;
        geomlib_spmat_from_triplets(rows, cols, vals, 4, 3, 2, &a);
        double b[3] = {1, 2, 3.5};
        double x[2] = {0, 0};
        check(geomlib_solve_cgls(&a, b, x, NULL) == GEOMLIB_SPARSE_SUCCESS, "CGLS rc");
        /* Normal equations: [2 1;1 2] x = [4.5; 5.5] => x = [1.3333, 2.0833]. */
        check_near(x[0] + x[1], 10.0 / 3.0, 1e-6, "CGLS sum");
        check_near(2*x[0] + x[1], 4.5, 1e-6, "CGLS normal eq 0");
        check_near(x[0] + 2*x[1], 5.5, 1e-6, "CGLS normal eq 1");

        double dense[6] = {1,0, 0,1, 1,1};
        double xqr[2] = {0, 0};
        check(geomlib_solve_qr_dense(dense, 3, 2, b, xqr) == GEOMLIB_SPARSE_SUCCESS, "QR rc");
        check_near(xqr[0], x[0], 1e-6, "QR agrees with CGLS 0");
        check_near(xqr[1], x[1], 1e-6, "QR agrees with CGLS 1");
        geomlib_spmat_free(&a);
    }

    /* Symmetric eigensolver. */
    {
        const double a[9] = {2, 1, 0, 1, 2, 0, 0, 0, 3};
        double ev[3], vec[9];
        check(geomlib_eigen_sym_dense(a, 3, ev, vec) == GEOMLIB_SPARSE_SUCCESS, "eigen rc");
        check_near(ev[0], 1.0, 1e-9, "eigenvalue 0");
        check_near(ev[1], 3.0, 1e-9, "eigenvalue 1");
        check_near(ev[2], 3.0, 1e-9, "eigenvalue 2");
        /* A v = lambda v for the first eigenvector. */
        double av0 = 2*vec[0] + 1*vec[1];
        double av1 = 1*vec[0] + 2*vec[1];
        double av2 = 3*vec[2];
        check_near(av0, ev[0]*vec[0], 1e-8, "eigenvector 0");
        check_near(av1, ev[0]*vec[1], 1e-8, "eigenvector 1");
        check_near(av2, ev[0]*vec[2], 1e-8, "eigenvector 2");
    }

    /* Geometry operators on a single triangle. */
    {
        geomlib_trimesh mesh;
        geomlib_trimesh_init(&mesh, 3, 1);
        geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0,0,0), NULL, NULL);
        geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(1,0,0), NULL, NULL);
        geomlib_trimesh_add_vertex(&mesh, geomlib_vec3_create(0,1,0), NULL, NULL);
        geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
        geomlib_spmat l, mass;
        check(geomlib_cotmatrix(&mesh, &l) == GEOMLIB_SPARSE_SUCCESS, "cotmatrix rc");
        check(geomlib_massmatrix(&mesh, &mass) == GEOMLIB_SPARSE_SUCCESS, "mass rc");
        check(l.nnz == 9, "cotmatrix dense 3x3");
        double ones[3] = {1, 1, 1};
        double rowsum[3];
        geomlib_spmat_mult_vec(&l, ones, rowsum);
        for (int i = 0; i < 3; ++i) check_near(rowsum[i], 0.0, 1e-9, "Laplacian row sum");
        double mass_sum;
        double mones[3];
        geomlib_spmat_mult_vec(&mass, ones, mones);
        mass_sum = mones[0] + mones[1] + mones[2];
        check_near(mass_sum, 0.5, 1e-9, "mass matrix total area");
        geomlib_spmat_free(&l);
        geomlib_spmat_free(&mass);
        geomlib_trimesh_free(&mesh);
    }

    if (failures != 0) {
        std::printf("%d sparse test(s) failed\n", failures);
        return 1;
    }
    std::printf("All sparse algebra tests passed\n");
    return 0;
}
