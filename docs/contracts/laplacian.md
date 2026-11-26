# Contract: cotangent Laplacian

**Status:** implemented (Wave 3).

## Source (paper only)

- Pinkall and Polthier, "Computing Discrete Minimal and Conformal Surfaces",
  Experimental Mathematics 1993.

## Approach

Assemble cotangent weights for each mesh edge from the two opposite corners,
emit symmetric CSR triplets, set the diagonal to the negative sum of off-diagonal
weights, and accumulate lumped vertex areas. The apply helper computes
`(L + scale*M) x`.

## API / constraints

`geomlib_build_laplacian(mesh, out)` and `geomlib_laplacian_apply(l, x, scale, y)`.
The stiffness matrix is symmetric with a constant nullspace. C11, wasm32 safe.

## Acceptance

`tests/test_laplacian.cpp`: a 5×5 planar grid gives total lumped area 16.0000
and maximum `L·1` of 0.000e+00. WASM smoke checks the nullspace.

## Lean target

Not yet formalized.
