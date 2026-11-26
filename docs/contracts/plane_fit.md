# Contract: 3D least-squares plane fitting

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard PCA plane fit: the plane normal is the smallest eigenvector of the
  centred covariance matrix.

## Approach

Compute the centroid and centred covariance of the points, extract the smallest
eigenvector by cyclic Jacobi rotations, and report the centroid as a point on
the plane plus the RMS distance of all points to it.

## API / constraints

`geomlib_fit_plane(points, n, out)` needs at least 3 non-collinear points. C11,
wasm32 safe.

## Acceptance

`tests/test_plane_fit.cpp`: a horizontal point set gives normal (0,0,1) and
RMS 0; a tilted plane z = 0.5x + 1 gives normal (−0.447, 0, 0.894) with RMS 0.
WASM smoke checks the horizontal case.

## Lean target

Not yet formalized.
