# Contract: principal component analysis

**Status:** implemented (Wave 3).

## Source (paper only)

- Pearson, "On Lines and Planes of Closest Fit to Systems of Points in Space",
  Philosophical Magazine 1901.

## Approach

Compute the mean and covariance of the point set, diagonalize the symmetric
3×3 covariance by cyclic Jacobi rotations, and return eigenvalues in descending
order with orthonormal eigenvectors.

## API / constraints

`geomlib_compute_pca(points, n, out)` requires at least 2 points. C11, wasm32
safe.

## Acceptance

`tests/test_pca.cpp`: an elongated cloud recovers the x axis as the first
principal direction with eigenvalue 76.939 versus 0.079 for the second.
WASM smoke checks the dominant eigenvalue.

## Lean target

Not yet formalized.
