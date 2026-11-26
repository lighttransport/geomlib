# Contract: sphere–AABB overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard closest-point test for a sphere and an axis-aligned bounding box;
  see Ericson, "Real-Time Collision Detection".

## Approach

Clamp the sphere centre to the box to obtain the closest box point, compute
the Euclidean distance from the centre to it, and compare with the radius. The
signed `distance` (centre-to-box minus radius) gives penetration depth when
negative.

## API / constraints

`geomlib_sphere_aabb_overlap(center, radius, bmin, bmax, out)` returns 1 on
overlap or tangency. C11, wasm32 safe.

## Acceptance

`tests/test_sphere_aabb.cpp`: a centre one unit from the face gives signed
distance +0.5 (miss) for radius 0.5 and −0.5 (hit) for radius 1.5; the closest
box point is the face projection; an interior centre overlaps. WASM smoke
checks the penetration depth.

## Lean target

Not yet formalized.
