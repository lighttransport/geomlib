# Contract: 2D ellipse fitting

**Status:** implemented (Wave 3).

## Source (paper only)

- PCA initialization followed by nonlinear least squares on the parametric
  ellipse (centre, semi-axes, orientation).

## Approach

Initialize from the mean and covariance eigenvectors, then minimize the sum of
squared normalized radial residuals by finite-difference gradient descent with
backtracking line search. Semi-axes are kept positive and ordered.

## API / constraints

`geomlib_fit_ellipse(points, n, out)` requires at least 5 points. C11, wasm32
safe.

## Acceptance

`tests/test_ellipse_fit.cpp`: a rotated ellipse with centre (1, −2), semi-axes
3 and 1.5 and angle 0.4 is recovered exactly with zero residual. WASM smoke
checks a 2×1 ellipse.

## Lean target

Not yet formalized.
