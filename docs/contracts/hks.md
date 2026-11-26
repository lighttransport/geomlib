# Contract: heat kernel signature

**Status:** implemented (Wave 3).

## Source (paper only)

- Sun, Ovsjanikov, Guibas, "A Concise and Provably Informative Multi-Scale
  Signature Based on Heat Diffusion", Computer Graphics Forum 2009.

## Approach

Compute the lowest Laplace-Beltrami eigenpairs with the existing manifold
harmonics solver, evaluate `sum_k exp(-lambda_k t) phi_k(i)^2` at each requested
time, and normalize each time slice by its maximum.

## API / constraints

`geomlib_compute_hks(mesh, eigen_count, times, ntimes, out)` with `out` of size
`vertex_count * ntimes` in row-major order. Times must be positive. The
descriptor depends on the quality of the eigenbasis; degenerate pole fans can
produce strongly peaked values. C11, wasm32 safe.

## Acceptance

`tests/test_hks.cpp`: a sphere at times 0.1 and 1.0 gives 234 finite positive
values in [0,1] that decay with time. WASM smoke checks one time slice.

## Lean target

Not yet formalized.
