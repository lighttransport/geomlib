# Contract: sparse direct solver

**Status:** implemented (Wave 3).

## Source (paper only)

- Cuthill & McKee, "Reducing the Bandwidth of Sparse Symmetric Matrices", 1969.
- Standard banded Cholesky factorization for SPD systems.

## Approach

Apply reverse Cuthill-McKee ordering to reduce bandwidth, copy the required
lower-triangular band into a compact banded array, factor with Cholesky, and
solve by forward/back substitution before undoing the permutation.

## API / constraints

`geomlib_solve_sparse_direct(a, b, x, bandwidth_out)` for square SPD sparse
matrices. C11, wasm32 safe. RCM currently visits neighbours in natural order.

## Acceptance

`tests/test_sparse_direct.cpp`: a 64x64 SPD tridiagonal system solves with
bandwidth 1 and residual 9.8e-26. WASM smoke solves a small system.

## Lean target

Not yet formalized.
