# Contract: mean-value cage deformation

**Status:** implemented (Wave 3).

## Source (paper only)

- Floater, "Mean Value Coordinates", Computer Aided Geometric Design 2003.

## Approach

For a query point q inside a planar cage, compute the angles subtended by
consecutive cage edges at q and set
`w_i = (tan(alpha_{i-1}/2) + tan(alpha_i/2)) / |c_i - q|`, normalized to sum
to one. Deform q by the weighted sum of the deformed cage vertices.

## API / constraints

`geomlib_deform_mean_value2d(cage, deformed, n, points, count, out)` and
`geomlib_mean_value_weights2d`. Cage vertices are handled exactly. The method
assumes a convex cage for positive weights; points outside may give negative
weights. C11, wasm32 safe.

## Acceptance

`tests/test_meancoords.cpp`: a square cage reproduces identity, translation
and uniform scale exactly, maps cage vertices exactly, and weights sum to one
and stay positive. WASM smoke checks a translation.

## Lean target

Not yet formalized.
