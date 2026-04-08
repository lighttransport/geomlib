/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Linear System Solver Implementation
 */

#include "geomlib/geomlib_allocator.h"
#include "geomlib/geomlib_c11.h"
#include "geomlib/geomlib_solver.h"

/* ============================================================================
 * Gaussian Elimination with Partial Pivoting (float)
 * ============================================================================ */

int geomlib_solve_gauss(size_t n, float *A, float *b) {
    if (n == 0 || A == NULL || b == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* Forward elimination with partial pivoting */
    for (size_t k = 0; k < n; k++) {
        /* Find pivot */
        size_t max_row = k;
        float max_val = fabsf(A[k * n + k]);
        
        for (size_t i = k + 1; i < n; i++) {
            float val = fabsf(A[i * n + k]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }
        
        /* Check for singular matrix */
        if (max_val < 1e-10f) {
            return GEOMLIB_SOLVER_ERROR_SINGULAR;
        }
        
        /* Swap rows if necessary */
        if (max_row != k) {
            for (size_t j = 0; j < n; j++) {
                float tmp = A[k * n + j];
                A[k * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            float tmp = b[k];
            b[k] = b[max_row];
            b[max_row] = tmp;
        }
        
        /* Eliminate column */
        for (size_t i = k + 1; i < n; i++) {
            float factor = A[i * n + k] / A[k * n + k];
            for (size_t j = k; j < n; j++) {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }
    
    /* Back substitution */
    for (size_t k = n; k > 0; k--) {
        size_t i = k - 1;
        float sum = b[i];
        for (size_t j = i + 1; j < n; j++) {
            sum -= A[i * n + j] * b[j];
        }
        b[i] = sum / A[i * n + i];
    }
    
    return GEOMLIB_SOLVER_SUCCESS;
}

/* ============================================================================
 * Gaussian Elimination with Partial Pivoting (double)
 * ============================================================================ */

int geomlib_solve_gauss_d(size_t n, double *A, double *b) {
    if (n == 0 || A == NULL || b == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* Forward elimination with partial pivoting */
    for (size_t k = 0; k < n; k++) {
        /* Find pivot */
        size_t max_row = k;
        double max_val = fabs(A[k * n + k]);
        
        for (size_t i = k + 1; i < n; i++) {
            double val = fabs(A[i * n + k]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }
        
        /* Check for singular matrix */
        if (max_val < 1e-15) {
            return GEOMLIB_SOLVER_ERROR_SINGULAR;
        }
        
        /* Swap rows if necessary */
        if (max_row != k) {
            for (size_t j = 0; j < n; j++) {
                double tmp = A[k * n + j];
                A[k * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            double tmp = b[k];
            b[k] = b[max_row];
            b[max_row] = tmp;
        }
        
        /* Eliminate column */
        for (size_t i = k + 1; i < n; i++) {
            double factor = A[i * n + k] / A[k * n + k];
            for (size_t j = k; j < n; j++) {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }
    
    /* Back substitution */
    for (size_t k = n; k > 0; k--) {
        size_t i = k - 1;
        double sum = b[i];
        for (size_t j = i + 1; j < n; j++) {
            sum -= A[i * n + j] * b[j];
        }
        b[i] = sum / A[i * n + i];
    }
    
    return GEOMLIB_SOLVER_SUCCESS;
}

/* ============================================================================
 * LU Decomposition with Partial Pivoting (float)
 * ============================================================================ */

int geomlib_lu_decompose(size_t n, float *A, int *pivot) {
    if (n == 0 || A == NULL || pivot == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* Initialize pivot array */
    for (size_t i = 0; i < n; i++) {
        pivot[i] = (int)i;
    }
    
    for (size_t k = 0; k < n; k++) {
        /* Find pivot */
        size_t max_row = k;
        float max_val = fabsf(A[k * n + k]);
        
        for (size_t i = k + 1; i < n; i++) {
            float val = fabsf(A[i * n + k]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }
        
        /* Check for singular matrix */
        if (max_val < 1e-10f) {
            return GEOMLIB_SOLVER_ERROR_SINGULAR;
        }
        
        /* Swap rows if necessary */
        if (max_row != k) {
            for (size_t j = 0; j < n; j++) {
                float tmp = A[k * n + j];
                A[k * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            int tmp = pivot[k];
            pivot[k] = pivot[max_row];
            pivot[max_row] = tmp;
        }
        
        /* Compute L and U */
        for (size_t i = k + 1; i < n; i++) {
            A[i * n + k] /= A[k * n + k];
            for (size_t j = k + 1; j < n; j++) {
                A[i * n + j] -= A[i * n + k] * A[k * n + j];
            }
        }
    }
    
    return GEOMLIB_SOLVER_SUCCESS;
}

int geomlib_lu_solve(size_t n, const float *LU, const int *pivot, float *b) {
    if (n == 0 || LU == NULL || pivot == NULL || b == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* Apply permutation to b */
    float *temp = (float *)geomlib_malloc(n * sizeof(float));
    if (temp == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    for (size_t i = 0; i < n; i++) {
        temp[i] = b[pivot[i]];
    }
    geomlib_memcpy(b, temp, n * sizeof(float));
    geomlib_free(temp);
    
    /* Forward substitution (solve Ly = Pb) */
    for (size_t i = 1; i < n; i++) {
        float sum = b[i];
        for (size_t j = 0; j < i; j++) {
            sum -= LU[i * n + j] * b[j];
        }
        b[i] = sum;
    }
    
    /* Back substitution (solve Ux = y) */
    for (size_t k = n; k > 0; k--) {
        size_t i = k - 1;
        float sum = b[i];
        for (size_t j = i + 1; j < n; j++) {
            sum -= LU[i * n + j] * b[j];
        }
        b[i] = sum / LU[i * n + i];
    }
    
    return GEOMLIB_SOLVER_SUCCESS;
}

/* ============================================================================
 * LU Decomposition with Partial Pivoting (double)
 * ============================================================================ */

int geomlib_lu_decompose_d(size_t n, double *A, int *pivot) {
    if (n == 0 || A == NULL || pivot == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* Initialize pivot array */
    for (size_t i = 0; i < n; i++) {
        pivot[i] = (int)i;
    }
    
    for (size_t k = 0; k < n; k++) {
        /* Find pivot */
        size_t max_row = k;
        double max_val = fabs(A[k * n + k]);
        
        for (size_t i = k + 1; i < n; i++) {
            double val = fabs(A[i * n + k]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }
        
        /* Check for singular matrix */
        if (max_val < 1e-15) {
            return GEOMLIB_SOLVER_ERROR_SINGULAR;
        }
        
        /* Swap rows if necessary */
        if (max_row != k) {
            for (size_t j = 0; j < n; j++) {
                double tmp = A[k * n + j];
                A[k * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            int tmp = pivot[k];
            pivot[k] = pivot[max_row];
            pivot[max_row] = tmp;
        }
        
        /* Compute L and U */
        for (size_t i = k + 1; i < n; i++) {
            A[i * n + k] /= A[k * n + k];
            for (size_t j = k + 1; j < n; j++) {
                A[i * n + j] -= A[i * n + k] * A[k * n + j];
            }
        }
    }
    
    return GEOMLIB_SOLVER_SUCCESS;
}

int geomlib_lu_solve_d(size_t n, const double *LU, const int *pivot, double *b) {
    if (n == 0 || LU == NULL || pivot == NULL || b == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* Apply permutation to b */
    double *temp = (double *)geomlib_malloc(n * sizeof(double));
    if (temp == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    for (size_t i = 0; i < n; i++) {
        temp[i] = b[pivot[i]];
    }
    geomlib_memcpy(b, temp, n * sizeof(double));
    geomlib_free(temp);
    
    /* Forward substitution (solve Ly = Pb) */
    for (size_t i = 1; i < n; i++) {
        double sum = b[i];
        for (size_t j = 0; j < i; j++) {
            sum -= LU[i * n + j] * b[j];
        }
        b[i] = sum;
    }
    
    /* Back substitution (solve Ux = y) */
    for (size_t k = n; k > 0; k--) {
        size_t i = k - 1;
        double sum = b[i];
        for (size_t j = i + 1; j < n; j++) {
            sum -= LU[i * n + j] * b[j];
        }
        b[i] = sum / LU[i * n + i];
    }
    
    return GEOMLIB_SOLVER_SUCCESS;
}

/* ============================================================================
 * Conjugate Gradient Solver (float)
 * ============================================================================ */

static void matvec_f(size_t n, const float *A, const float *x, float *y) {
    for (size_t i = 0; i < n; i++) {
        float sum = 0.0f;
        for (size_t j = 0; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }
        y[i] = sum;
    }
}

static float dot_f(size_t n, const float *a, const float *b) {
    float sum = 0.0f;
    for (size_t i = 0; i < n; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

int geomlib_solve_cg(size_t n, const float *A, const float *b, float *x,
                     int max_iter, float tol, int *iter_out) {
    if (n == 0 || A == NULL || b == NULL || x == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    float *r = (float *)geomlib_malloc(n * sizeof(float));
    float *p = (float *)geomlib_malloc(n * sizeof(float));
    float *Ap = (float *)geomlib_malloc(n * sizeof(float));
    
    if (r == NULL || p == NULL || Ap == NULL) {
        geomlib_free(r);
        geomlib_free(p);
        geomlib_free(Ap);
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* r = b - A*x */
    matvec_f(n, A, x, r);
    for (size_t i = 0; i < n; i++) {
        r[i] = b[i] - r[i];
        p[i] = r[i];
    }
    
    float rsold = dot_f(n, r, r);
    float tol_sq = tol * tol;
    
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        matvec_f(n, A, p, Ap);
        float pAp = dot_f(n, p, Ap);
        
        if (fabsf(pAp) < 1e-20f) {
            break;
        }
        
        float alpha = rsold / pAp;
        
        for (size_t i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        
        float rsnew = dot_f(n, r, r);
        
        if (rsnew < tol_sq) {
            iter++;
            break;
        }
        
        float beta = rsnew / rsold;
        for (size_t i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i];
        }
        
        rsold = rsnew;
    }
    
    if (iter_out != NULL) {
        *iter_out = iter;
    }
    
    geomlib_free(r);
    geomlib_free(p);
    geomlib_free(Ap);
    
    return (iter < max_iter) ? GEOMLIB_SOLVER_SUCCESS : GEOMLIB_SOLVER_ERROR_NO_CONVERGE;
}

/* ============================================================================
 * Conjugate Gradient Solver (double)
 * ============================================================================ */

static void matvec_d(size_t n, const double *A, const double *x, double *y) {
    for (size_t i = 0; i < n; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }
        y[i] = sum;
    }
}

static double dot_d(size_t n, const double *a, const double *b) {
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

int geomlib_solve_cg_d(size_t n, const double *A, const double *b, double *x,
                       int max_iter, double tol, int *iter_out) {
    if (n == 0 || A == NULL || b == NULL || x == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    double *r = (double *)geomlib_malloc(n * sizeof(double));
    double *p = (double *)geomlib_malloc(n * sizeof(double));
    double *Ap = (double *)geomlib_malloc(n * sizeof(double));
    
    if (r == NULL || p == NULL || Ap == NULL) {
        geomlib_free(r);
        geomlib_free(p);
        geomlib_free(Ap);
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* r = b - A*x */
    matvec_d(n, A, x, r);
    for (size_t i = 0; i < n; i++) {
        r[i] = b[i] - r[i];
        p[i] = r[i];
    }
    
    double rsold = dot_d(n, r, r);
    double tol_sq = tol * tol;
    
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        matvec_d(n, A, p, Ap);
        double pAp = dot_d(n, p, Ap);
        
        if (fabs(pAp) < 1e-30) {
            break;
        }
        
        double alpha = rsold / pAp;
        
        for (size_t i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        
        double rsnew = dot_d(n, r, r);
        
        if (rsnew < tol_sq) {
            iter++;
            break;
        }
        
        double beta = rsnew / rsold;
        for (size_t i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i];
        }
        
        rsold = rsnew;
    }
    
    if (iter_out != NULL) {
        *iter_out = iter;
    }
    
    geomlib_free(r);
    geomlib_free(p);
    geomlib_free(Ap);
    
    return (iter < max_iter) ? GEOMLIB_SOLVER_SUCCESS : GEOMLIB_SOLVER_ERROR_NO_CONVERGE;
}

/* ============================================================================
 * Jacobi Iterative Solver (float)
 * ============================================================================ */

int geomlib_solve_jacobi(size_t n, const float *A, const float *b, float *x,
                         int max_iter, float tol, int *iter_out) {
    if (n == 0 || A == NULL || b == NULL || x == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    float *x_new = (float *)geomlib_malloc(n * sizeof(float));
    if (x_new == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        float max_diff = 0.0f;
        
        for (size_t i = 0; i < n; i++) {
            if (fabsf(A[i * n + i]) < 1e-20f) {
                geomlib_free(x_new);
                return GEOMLIB_SOLVER_ERROR_SINGULAR;
            }
            
            float sum = b[i];
            for (size_t j = 0; j < n; j++) {
                if (j != i) {
                    sum -= A[i * n + j] * x[j];
                }
            }
            x_new[i] = sum / A[i * n + i];
            
            float diff = fabsf(x_new[i] - x[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
        
        geomlib_memcpy(x, x_new, n * sizeof(float));
        
        if (max_diff < tol) {
            iter++;
            break;
        }
    }
    
    if (iter_out != NULL) {
        *iter_out = iter;
    }
    
    geomlib_free(x_new);
    
    return (iter < max_iter) ? GEOMLIB_SOLVER_SUCCESS : GEOMLIB_SOLVER_ERROR_NO_CONVERGE;
}

/* ============================================================================
 * Jacobi Iterative Solver (double)
 * ============================================================================ */

int geomlib_solve_jacobi_d(size_t n, const double *A, const double *b, double *x,
                           int max_iter, double tol, int *iter_out) {
    if (n == 0 || A == NULL || b == NULL || x == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    double *x_new = (double *)geomlib_malloc(n * sizeof(double));
    if (x_new == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        double max_diff = 0.0;
        
        for (size_t i = 0; i < n; i++) {
            if (fabs(A[i * n + i]) < 1e-30) {
                geomlib_free(x_new);
                return GEOMLIB_SOLVER_ERROR_SINGULAR;
            }
            
            double sum = b[i];
            for (size_t j = 0; j < n; j++) {
                if (j != i) {
                    sum -= A[i * n + j] * x[j];
                }
            }
            x_new[i] = sum / A[i * n + i];
            
            double diff = fabs(x_new[i] - x[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
        
        geomlib_memcpy(x, x_new, n * sizeof(double));
        
        if (max_diff < tol) {
            iter++;
            break;
        }
    }
    
    if (iter_out != NULL) {
        *iter_out = iter;
    }
    
    geomlib_free(x_new);
    
    return (iter < max_iter) ? GEOMLIB_SOLVER_SUCCESS : GEOMLIB_SOLVER_ERROR_NO_CONVERGE;
}

/* ============================================================================
 * Power Iteration for Eigenvalues (float)
 * ============================================================================ */

int geomlib_eigen_power(size_t n, const float *A, float *eigenvalue, float *eigenvector,
                        int max_iter, float tol) {
    if (n == 0 || A == NULL || eigenvalue == NULL || eigenvector == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    float *y = (float *)geomlib_malloc(n * sizeof(float));
    if (y == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* Initialize eigenvector */
    for (size_t i = 0; i < n; i++) {
        eigenvector[i] = 1.0f;
    }
    
    float lambda_old = 0.0f;
    
    for (int iter = 0; iter < max_iter; iter++) {
        /* y = A * x */
        matvec_f(n, A, eigenvector, y);
        
        /* Find max component for normalization */
        float max_val = 0.0f;
        for (size_t i = 0; i < n; i++) {
            if (fabsf(y[i]) > max_val) {
                max_val = fabsf(y[i]);
            }
        }
        
        if (max_val < 1e-20f) {
            geomlib_free(y);
            return GEOMLIB_SOLVER_ERROR_SINGULAR;
        }
        
        /* Normalize */
        float lambda = max_val;
        for (size_t i = 0; i < n; i++) {
            eigenvector[i] = y[i] / max_val;
        }
        
        /* Check convergence */
        if (fabsf(lambda - lambda_old) < tol) {
            *eigenvalue = lambda;
            geomlib_free(y);
            return GEOMLIB_SOLVER_SUCCESS;
        }
        
        lambda_old = lambda;
    }
    
    *eigenvalue = lambda_old;
    geomlib_free(y);
    return GEOMLIB_SOLVER_ERROR_NO_CONVERGE;
}

/* ============================================================================
 * Power Iteration for Eigenvalues (double)
 * ============================================================================ */

int geomlib_eigen_power_d(size_t n, const double *A, double *eigenvalue, double *eigenvector,
                          int max_iter, double tol) {
    if (n == 0 || A == NULL || eigenvalue == NULL || eigenvector == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    double *y = (double *)geomlib_malloc(n * sizeof(double));
    if (y == NULL) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* Initialize eigenvector */
    for (size_t i = 0; i < n; i++) {
        eigenvector[i] = 1.0;
    }
    
    double lambda_old = 0.0;
    
    for (int iter = 0; iter < max_iter; iter++) {
        /* y = A * x */
        matvec_d(n, A, eigenvector, y);
        
        /* Find max component for normalization */
        double max_val = 0.0;
        for (size_t i = 0; i < n; i++) {
            if (fabs(y[i]) > max_val) {
                max_val = fabs(y[i]);
            }
        }
        
        if (max_val < 1e-30) {
            geomlib_free(y);
            return GEOMLIB_SOLVER_ERROR_SINGULAR;
        }
        
        /* Normalize */
        double lambda = max_val;
        for (size_t i = 0; i < n; i++) {
            eigenvector[i] = y[i] / max_val;
        }
        
        /* Check convergence */
        if (fabs(lambda - lambda_old) < tol) {
            *eigenvalue = lambda;
            geomlib_free(y);
            return GEOMLIB_SOLVER_SUCCESS;
        }
        
        lambda_old = lambda;
    }
    
    *eigenvalue = lambda_old;
    geomlib_free(y);
    return GEOMLIB_SOLVER_ERROR_NO_CONVERGE;
}

/* ============================================================================
 * Least Squares Solver (float)
 * ============================================================================ */

int geomlib_solve_least_squares(size_t m, size_t n, float *A, float *b, float *x) {
    if (m == 0 || n == 0 || A == NULL || b == NULL || x == NULL || m < n) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* Compute A^T * A */
    float *ATA = (float *)geomlib_malloc(n * n * sizeof(float));
    /* Compute A^T * b */
    float *ATb = (float *)geomlib_malloc(n * sizeof(float));
    
    if (ATA == NULL || ATb == NULL) {
        geomlib_free(ATA);
        geomlib_free(ATb);
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* ATA = A^T * A */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            float sum = 0.0f;
            for (size_t k = 0; k < m; k++) {
                sum += A[k * n + i] * A[k * n + j];
            }
            ATA[i * n + j] = sum;
        }
    }
    
    /* ATb = A^T * b */
    for (size_t i = 0; i < n; i++) {
        float sum = 0.0f;
        for (size_t k = 0; k < m; k++) {
            sum += A[k * n + i] * b[k];
        }
        ATb[i] = sum;
    }
    
    /* Solve ATA * x = ATb using Gaussian elimination */
    int result = geomlib_solve_gauss(n, ATA, ATb);
    
    if (result == GEOMLIB_SOLVER_SUCCESS) {
        geomlib_memcpy(x, ATb, n * sizeof(float));
    }
    
    geomlib_free(ATA);
    geomlib_free(ATb);
    
    return result;
}

/* ============================================================================
 * Least Squares Solver (double)
 * ============================================================================ */

int geomlib_solve_least_squares_d(size_t m, size_t n, double *A, double *b, double *x) {
    if (m == 0 || n == 0 || A == NULL || b == NULL || x == NULL || m < n) {
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* Compute A^T * A */
    double *ATA = (double *)geomlib_malloc(n * n * sizeof(double));
    /* Compute A^T * b */
    double *ATb = (double *)geomlib_malloc(n * sizeof(double));
    
    if (ATA == NULL || ATb == NULL) {
        geomlib_free(ATA);
        geomlib_free(ATb);
        return GEOMLIB_SOLVER_ERROR_INVALID;
    }
    
    /* ATA = A^T * A */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < m; k++) {
                sum += A[k * n + i] * A[k * n + j];
            }
            ATA[i * n + j] = sum;
        }
    }
    
    /* ATb = A^T * b */
    for (size_t i = 0; i < n; i++) {
        double sum = 0.0;
        for (size_t k = 0; k < m; k++) {
            sum += A[k * n + i] * b[k];
        }
        ATb[i] = sum;
    }
    
    /* Solve ATA * x = ATb using Gaussian elimination */
    int result = geomlib_solve_gauss_d(n, ATA, ATb);
    
    if (result == GEOMLIB_SOLVER_SUCCESS) {
        geomlib_memcpy(x, ATb, n * sizeof(double));
    }
    
    geomlib_free(ATA);
    geomlib_free(ATb);
    
    return result;
}
