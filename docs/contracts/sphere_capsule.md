# Contract: sphere–capsule overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Sphere-capsule collision as the closest-point distance between the sphere
  centre and the capsule axis compared with the sum of radii.

## Approach

Project the sphere centre onto the capsule axis with the existing capsule point
distance helper (radius zero), then subtract the sphere and capsule radii. The
residue is negative when the sphere penetrates the capsule.

## API / constraints

`geomlib_sphere_capsule_overlap(center, radius, a, b, cap_radius, out)` requires
non-negative radii. C11, wasm32 safe.

## Acceptance

`tests/test_sphere_capsule.cpp`: a centre one unit from the axis with radii
summing to one is tangent; a smaller sphere misses by 0.3; a centre on the axis
penetrates by 1. WASM smoke checks the penetrating case.

## Lean target

Not yet formalized.
