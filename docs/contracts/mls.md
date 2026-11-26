# Contract: MLS surface projection

**Status:** implemented (Wave 3).

## Source (paper only)

- Alexa, Behr, Cohen-Or, Fleishman, Silva, "Computing and Rendering Point Set
  Surfaces", IEEE TVCG 2003 (moving least squares projection).

## Approach

Find k nearest neighbours with a k-d tree, fit a Gaussian-weighted local plane
by PCA (cyclic Jacobi eigen decomposition of the 3×3 covariance), project the
query onto that plane, and repeat for the requested number of iterations. The
smallest eigenvector gives the local surface normal.

## API / constraints

`geomlib_mls_project(cloud, query, params, projected, normal)` with
`geomlib_mls_params_default`. A Gaussian width of 0 uses the farthest neighbour
distance. C11, wasm32 safe.

## Acceptance

`tests/test_mls.cpp`: a planar cloud projects a query from z=0.3 to z=0.0000
with vertical normal; a spherical cloud projects a radial query from distance
0.2130 to 0.0799 from the nearest sample, with radius 0.96 and unit normal.
WASM smoke checks the planar case.

## Lean target

Not yet formalized.
