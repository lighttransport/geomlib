# Contract: capsule collision

**Status:** implemented (Wave 3).

## Source (paper only)

- Capsule as the Minkowski sum of a segment and a sphere; standard closest-point
  collision test.

## Approach

Delegate to the existing finite segment closest-point solver and compare the
axis distance with the sum of radii; the signed residue is the penetration
depth. Point-to-capsule distance clamps the query onto the axis and subtracts
the radius.

## API / constraints

`geomlib_capsule_capsule_overlap(a0,a1,r0,b0,b1,r1,out)` and
`geomlib_capsule_point_distance(q,a0,a1,r,axis_point)` require non-negative
radii. C11, wasm32 safe.

## Acceptance

`tests/test_capsule.cpp`: collinear capsules with a one-unit gap are tangent at
distance 0; a larger gap reports 2; perpendicular crossing axes penetrate by 1;
a point two units out gives surface distance 1.5. WASM smoke checks tangency.

## Note

The shared `geomlib_segment_closest_points` parallel branch was corrected to
try all endpoint-to-segment projections; its original parallel test still
passes and collinear capsules now work.

## Lean target

Not yet formalized.
