# Contract: rigid point-set fitting

**Status:** implemented (Wave 3).

## Source (paper only)

- Horn, "Closed-Form Solution of Absolute Orientation Using Unit Quaternions",
  JOSA 1987.

## Approach

Compute source/target centroids, build the cross-covariance matrix, assemble
Horn's symmetric 4×4 N matrix, and extract the largest eigenvector by cyclic
Jacobi rotations to get the rotation quaternion. Translation follows from the
centroids, and optional uniform scale from the projected covariance ratio.

## API / constraints

`geomlib_fit_rigid_points(src, dst, count, allow_scale, transform, rms_out)`
requires at least 3 correspondences. `transform` is column-major. C11, wasm32
safe.

## Acceptance

`tests/test_rigid_fit.cpp`: a known rotation + translation is recovered with
RMS 0.000001; enabling scale recovers a 2× similarity with the same RMS. WASM
smoke checks a translation fit.

## Lean target

Not yet formalized.
