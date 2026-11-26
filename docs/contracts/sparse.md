# Contract: sparse linear algebra

**Status:** implemented (Wave 0). **Sources:** standard numerical methods —
conjugate gradients (Hestenes–Stiefel), CGLS (Paige–Saunders), Householder QR,
cyclic Jacobi eigensolver. No external library source consulted.

## Specification

- CSR construction from triplets with duplicate coalescing and sorted columns.
- Jacobi-preconditioned CG for square SPD systems.
- CGLS for min ||Ax−b||; convergence measured on ||Aᵀr||, not ||b−Ax||.
- Dense Householder QR least squares and symmetric Jacobi eigensolver.
- `geomlib_cotmatrix` / `geomlib_massmatrix`: double-precision discrete
  Laplace–Beltrami stiffness and lumped area mass matrices.

## Derivation

Formulas from standard numerical analysis texts; no code copied. The cotangent
weights use the standard cot(angle at opposite vertex)/2 convention.

## API / constraints

C11, `double`, malloc-backed CSR arrays with explicit `geomlib_spmat_free`.
WASM: yes (32-bit `size_t` is sufficient for interactive meshes).

## Acceptance

`tests/test_sparse.cpp`: duplicate coalescing, matvec, CG solve, CGLS normal
equations, QR/CGLS agreement, symmetric eigensystem, Laplacian row-sum zero,
mass matrix total area.

## Lean target

Discrete operator identities (`d∘d=0`, cotan Laplacian = div∘grad) planned for
the full formalization layer; not yet proved.
