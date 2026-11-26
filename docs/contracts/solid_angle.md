# Contract: solid angle and winding number

**Status:** implemented (Wave 3).

## Source (paper only)

- Van Oosterom and Strackee, "The Solid Angle of a Plane Triangle", IEEE
  Transactions on Biomedical Engineering 1983.

## Approach

For vectors a, b, c from the query point to the triangle vertices, compute
`Omega = 2 atan2(a . (b x c), |a||b||c| + (a.b)|c| + (b.c)|a| + (c.a)|b|)`,
which is the signed solid angle. Polygon fans and closed-mesh sums accumulate
it; winding number is the sum divided by 4pi.

## API / constraints

`geomlib_solid_angle_triangle`, `geomlib_solid_angle_polygon` and
`geomlib_solid_angle_mesh` return steradians. `geomlib_winding_number_mesh`
returns the generalised winding number. C11, wasm32 safe.

## Acceptance

`tests/test_solid_angle.cpp`: an octant triangle subtends pi/2 exactly,
reversing winding flips the sign, a closed box gives 4pi inside (winding 1) and
0 outside. WASM smoke checks the octant value and interior winding.

## Lean target

Not yet formalized.
