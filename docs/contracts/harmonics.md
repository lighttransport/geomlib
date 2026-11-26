# Contract: manifold harmonics

**Status:** implemented (Wave 3).

## Sources (papers only)

- Desbrun, Meyer, Schroder, Barr, "Implicit fairing of irregular meshes using
  diffusion and curvature flow", SIGGRAPH 1999.
- Levy, "Laplace-Beltrami eigenfunctions towards an algorithm that understands
  geometry", SGP 2006.

## Approach

Form the generalized symmetric eigenproblem L x = lambda M x with cotangent
stiffness and lumped area mass; reduce it to the standard symmetric problem
M^{-1/2} L M^{-1/2} y = lambda y and solve with the existing cyclic Jacobi
eigensolver. Eigenvectors are transformed back with M^{-1/2}. Vertices with
zero lumped area (isolated/degenerate) receive a small regularized mass.

## API / constraints

`geomlib_compute_manifold_harmonics(mesh, count, out)`. Supports meshes up to
512 vertices (dense eigensolver). C11, wasm32 safe.

## Acceptance

`tests/test_harmonics.cpp`: sphere spectrum starts at zero, is ascending, and
the first two modes are mass-orthogonal. WASM smoke checks the zero eigenvalue.

## Lean target

Spectral convergence is future work.
