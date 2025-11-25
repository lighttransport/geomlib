/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Linear System Solvers
 * Core functions implemented in C for maximum portability
 */

#ifndef GEOMLIB_SOLVER_H
#define GEOMLIB_SOLVER_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Return Codes
 * ============================================================================ */

#define GEOMLIB_SOLVER_SUCCESS           0
#define GEOMLIB_SOLVER_ERROR_SINGULAR    1
#define GEOMLIB_SOLVER_ERROR_NO_CONVERGE 2
#define GEOMLIB_SOLVER_ERROR_INVALID     3

/* ============================================================================
 * Dense Linear System Solvers
 * Solve Ax = b where A is n x n matrix, x and b are n-vectors
 * ============================================================================ */

/**
 * Solve linear system using Gaussian elimination with partial pivoting
 * 
 * @param n     Size of the system (n x n matrix)
 * @param A     Input matrix (row-major, will be modified)
 * @param b     Right-hand side vector (will contain solution on success)
 * @return      GEOMLIB_SOLVER_SUCCESS on success, error code otherwise
 */
int geomlib_solve_gauss(size_t n, float *A, float *b);
int geomlib_solve_gauss_d(size_t n, double *A, double *b);

/**
 * LU decomposition with partial pivoting
 * A = P * L * U where P is permutation matrix
 * 
 * @param n     Size of the matrix
 * @param A     Input matrix (row-major, will contain L and U on output)
 * @param pivot Permutation array (size n)
 * @return      GEOMLIB_SOLVER_SUCCESS on success, error code otherwise
 */
int geomlib_lu_decompose(size_t n, float *A, int *pivot);
int geomlib_lu_decompose_d(size_t n, double *A, int *pivot);

/**
 * Solve using existing LU decomposition
 * 
 * @param n     Size of the system
 * @param LU    LU decomposition from geomlib_lu_decompose
 * @param pivot Permutation array from geomlib_lu_decompose
 * @param b     Right-hand side (will contain solution on success)
 * @return      GEOMLIB_SOLVER_SUCCESS on success
 */
int geomlib_lu_solve(size_t n, const float *LU, const int *pivot, float *b);
int geomlib_lu_solve_d(size_t n, const double *LU, const int *pivot, double *b);

/* ============================================================================
 * Iterative Solvers (for sparse/large systems)
 * ============================================================================ */

/**
 * Conjugate Gradient solver for symmetric positive-definite systems
 * 
 * @param n         Size of the system
 * @param A         Matrix (row-major, n x n)
 * @param b         Right-hand side vector
 * @param x         Initial guess (input) and solution (output)
 * @param max_iter  Maximum iterations
 * @param tol       Convergence tolerance
 * @param iter_out  Actual iterations performed (can be NULL)
 * @return          GEOMLIB_SOLVER_SUCCESS on success
 */
int geomlib_solve_cg(size_t n, const float *A, const float *b, float *x,
                     int max_iter, float tol, int *iter_out);
int geomlib_solve_cg_d(size_t n, const double *A, const double *b, double *x,
                       int max_iter, double tol, int *iter_out);

/**
 * Jacobi iterative solver
 * 
 * @param n         Size of the system
 * @param A         Matrix (row-major, n x n)
 * @param b         Right-hand side vector
 * @param x         Initial guess (input) and solution (output)
 * @param max_iter  Maximum iterations
 * @param tol       Convergence tolerance
 * @param iter_out  Actual iterations performed (can be NULL)
 * @return          GEOMLIB_SOLVER_SUCCESS on success
 */
int geomlib_solve_jacobi(size_t n, const float *A, const float *b, float *x,
                         int max_iter, float tol, int *iter_out);
int geomlib_solve_jacobi_d(size_t n, const double *A, const double *b, double *x,
                           int max_iter, double tol, int *iter_out);

/* ============================================================================
 * Eigenvalue/Eigenvector Computation
 * ============================================================================ */

/**
 * Power iteration to find dominant eigenvalue and eigenvector
 * 
 * @param n            Size of the matrix
 * @param A            Input matrix (row-major, n x n)
 * @param eigenvalue   Output eigenvalue
 * @param eigenvector  Output eigenvector (size n)
 * @param max_iter     Maximum iterations
 * @param tol          Convergence tolerance
 * @return             GEOMLIB_SOLVER_SUCCESS on success
 */
int geomlib_eigen_power(size_t n, const float *A, float *eigenvalue, float *eigenvector,
                        int max_iter, float tol);
int geomlib_eigen_power_d(size_t n, const double *A, double *eigenvalue, double *eigenvector,
                          int max_iter, double tol);

/* ============================================================================
 * Least Squares Solver
 * ============================================================================ */

/**
 * Solve overdetermined system Ax = b in least squares sense
 * Uses normal equations: A^T A x = A^T b
 * 
 * @param m     Number of rows in A (equations)
 * @param n     Number of columns in A (unknowns)
 * @param A     Matrix (row-major, m x n, will be modified)
 * @param b     Right-hand side (size m, will be modified)
 * @param x     Solution vector (size n)
 * @return      GEOMLIB_SOLVER_SUCCESS on success
 */
int geomlib_solve_least_squares(size_t m, size_t n, float *A, float *b, float *x);
int geomlib_solve_least_squares_d(size_t m, size_t n, double *A, double *b, double *x);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SOLVER_H */
