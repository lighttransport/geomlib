# Contract: plane–sphere intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Analytic intersection of a plane and a sphere in 3D Euclidean space.

## Approach

Normalize the plane normal, compute the signed distance from the sphere centre
to the plane, and derive the circle radius by Pythagoras with the circle centre
at the perpendicular projection of the sphere centre.

## API / constraints

`geomlib_plane_sphere_intersect(plane_point, plane_normal, sphere_center,
radius)` returns none / tangent / circle. C11, wasm32 safe.

## Acceptance

`tests/test_plane_sphere.cpp`: a centred sphere of radius 2 gives a great circle
of radius 2; a centre one unit above gives √3; a centre two units above is
tangent; three units above is disjoint. WASM smoke checks the √3 circle.

## Lean target

Not yet formalized.
