# Contract: sphere–triangle overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Closest-point-on-triangle collision test; see Ericson, "Real-Time Collision
  Detection".

## Approach

Project the sphere centre onto the triangle with the existing closest-point
helper, then compare the centre-to-triangle distance with the radius. The signed
residue reports penetration depth.

## API / constraints

`geomlib_sphere_triangle_overlap(center, radius, v0, v1, v2, out)` requires a
non-negative radius. C11, wasm32 safe.

## Acceptance

`tests/test_sphere_triangle.cpp`: a centre 0.1 above the interior with radius
0.2 gives penetration 0.1 at the projected interior point; a centre near the
hypotenuse edge gives penetration 0.5 with the closest point on x+y=2; a distant
sphere misses. WASM smoke checks the interior penetration.

## Lean target

Not yet formalized.
