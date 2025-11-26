/*
 * geomlib - Solver Tests
 */

#include <cstdio>
#include <cmath>
#include <cstring>
#include "geomlib/geomlib_solver.h"

static int test_count = 0;
static int pass_count = 0;

#define EPSILON 1e-4f

static bool float_eq(float a, float b, float eps = EPSILON) {
    return std::fabs(a - b) < eps;
}

static void test(bool condition, const char* name) {
    test_count++;
    if (condition) {
        pass_count++;
        std::printf("[PASS] %s\n", name);
    } else {
        std::printf("[FAIL] %s\n", name);
    }
}

static void test_gauss_2x2() {
    // Solve: 2x + y = 5
    //        x + 3y = 7
    // Solution: x = 1.6, y = 1.8
    
    float A[] = {2.0f, 1.0f,
                 1.0f, 3.0f};
    float b[] = {5.0f, 7.0f};
    
    int result = geomlib_solve_gauss(2, A, b);
    
    test(result == GEOMLIB_SOLVER_SUCCESS, "gauss_2x2_success");
    test(float_eq(b[0], 1.6f) && float_eq(b[1], 1.8f), "gauss_2x2_solution");
}

static void test_gauss_3x3() {
    // Solve: x + y + z = 6
    //        2x + y - z = 1  
    //        3x - y + z = 4
    // Solution: x = 1, y = 2, z = 3
    
    float A[] = {1.0f, 1.0f, 1.0f,
                 2.0f, 1.0f, -1.0f,
                 3.0f, -1.0f, 1.0f};
    float b[] = {6.0f, 1.0f, 4.0f};
    
    int result = geomlib_solve_gauss(3, A, b);
    
    test(result == GEOMLIB_SOLVER_SUCCESS, "gauss_3x3_success");
    test(float_eq(b[0], 1.0f) && float_eq(b[1], 2.0f) && float_eq(b[2], 3.0f), "gauss_3x3_solution");
}

static void test_gauss_singular() {
    // Singular matrix
    float A[] = {1.0f, 2.0f,
                 2.0f, 4.0f};  // Second row is 2x first row
    float b[] = {3.0f, 6.0f};
    
    int result = geomlib_solve_gauss(2, A, b);
    test(result == GEOMLIB_SOLVER_ERROR_SINGULAR, "gauss_singular");
}

static void test_lu_decompose() {
    float A[] = {2.0f, 1.0f,
                 1.0f, 3.0f};
    int pivot[2];
    
    int result = geomlib_lu_decompose(2, A, pivot);
    test(result == GEOMLIB_SOLVER_SUCCESS, "lu_decompose_success");
}

static void test_lu_solve() {
    // Same system as gauss_2x2
    float A[] = {2.0f, 1.0f,
                 1.0f, 3.0f};
    float b[] = {5.0f, 7.0f};
    int pivot[2];
    
    int result = geomlib_lu_decompose(2, A, pivot);
    test(result == GEOMLIB_SOLVER_SUCCESS, "lu_solve_decompose");
    
    result = geomlib_lu_solve(2, A, pivot, b);
    test(result == GEOMLIB_SOLVER_SUCCESS, "lu_solve_success");
    test(float_eq(b[0], 1.6f) && float_eq(b[1], 1.8f), "lu_solve_solution");
}

static void test_cg_spd() {
    // Symmetric positive definite matrix
    // A = [4 1; 1 3]
    float A[] = {4.0f, 1.0f,
                 1.0f, 3.0f};
    float b[] = {1.0f, 2.0f};
    float x[] = {0.0f, 0.0f};
    int iter;
    
    int result = geomlib_solve_cg(2, A, b, x, 100, 1e-6f, &iter);
    
    // Verify: A*x = b
    float Ax0 = A[0] * x[0] + A[1] * x[1];
    float Ax1 = A[2] * x[0] + A[3] * x[1];
    
    test(result == GEOMLIB_SOLVER_SUCCESS, "cg_spd_success");
    test(float_eq(Ax0, b[0], 1e-3f) && float_eq(Ax1, b[1], 1e-3f), "cg_spd_solution");
}

static void test_jacobi() {
    // Diagonally dominant matrix
    float A[] = {4.0f, 1.0f,
                 1.0f, 3.0f};
    float b[] = {1.0f, 2.0f};
    float x[] = {0.0f, 0.0f};
    int iter;
    
    int result = geomlib_solve_jacobi(2, A, b, x, 100, 1e-6f, &iter);
    
    // Verify solution
    float Ax0 = A[0] * x[0] + A[1] * x[1];
    float Ax1 = A[2] * x[0] + A[3] * x[1];
    
    test(result == GEOMLIB_SOLVER_SUCCESS, "jacobi_success");
    test(float_eq(Ax0, b[0], 1e-3f) && float_eq(Ax1, b[1], 1e-3f), "jacobi_solution");
}

static void test_power_iteration() {
    // Matrix with known dominant eigenvalue
    // A = [2 0; 0 1] has eigenvalue 2 with eigenvector [1, 0]
    float A[] = {2.0f, 0.0f,
                 0.0f, 1.0f};
    float eigenvalue;
    float eigenvector[2];
    
    int result = geomlib_eigen_power(2, A, &eigenvalue, eigenvector, 100, 1e-6f);
    
    test(result == GEOMLIB_SOLVER_SUCCESS, "power_iteration_success");
    test(float_eq(eigenvalue, 2.0f, 1e-3f), "power_iteration_eigenvalue");
}

static void test_least_squares() {
    // Overdetermined system (3 equations, 2 unknowns)
    // Fit a line y = mx + c to points (0,1), (1,2), (2,4)
    // A = [0 1; 1 1; 2 1], b = [1; 2; 4]
    
    float A[] = {0.0f, 1.0f,
                 1.0f, 1.0f,
                 2.0f, 1.0f};
    float b[] = {1.0f, 2.0f, 4.0f};
    float x[2];
    
    int result = geomlib_solve_least_squares(3, 2, A, b, x);
    
    test(result == GEOMLIB_SOLVER_SUCCESS, "least_squares_success");
    // The least squares solution should be approximately m=1.5, c=0.67
    test(float_eq(x[0], 1.5f, 0.1f), "least_squares_slope");
}

static void test_gauss_double() {
    // Double precision test
    double A[] = {2.0, 1.0,
                  1.0, 3.0};
    double b[] = {5.0, 7.0};
    
    int result = geomlib_solve_gauss_d(2, A, b);
    
    test(result == GEOMLIB_SOLVER_SUCCESS, "gauss_double_success");
    test(std::fabs(b[0] - 1.6) < 1e-10 && std::fabs(b[1] - 1.8) < 1e-10, "gauss_double_solution");
}

int main() {
    std::printf("=== Solver Tests ===\n\n");
    
    test_gauss_2x2();
    test_gauss_3x3();
    test_gauss_singular();
    test_lu_decompose();
    test_lu_solve();
    test_cg_spd();
    test_jacobi();
    test_power_iteration();
    test_least_squares();
    test_gauss_double();
    
    std::printf("\n=== Results: %d/%d tests passed ===\n", pass_count, test_count);
    
    return (pass_count == test_count) ? 0 : 1;
}
