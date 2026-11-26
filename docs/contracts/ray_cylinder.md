# Contract: ray–cylinder intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Analytic intersection of a ray with an infinite cylinder in 3D.

## Approach

Project the ray origin and direction onto the plane perpendicular to the axis,
solve the resulting quadratic for the cylinder radius, and return the nearest
non-negative root with its surface point and outward radial normal.

## API / constraints

`geomlib_ray_cylinder_intersect(origin, direction, axis_point, axis_dir, radius,
out)` returns 1 on a hit. Rays parallel to the axis are not reported. C11,
wasm32 safe.

## Acceptance

`tests/test_ray_cylinder.cpp`: a unit z-axis cylinder hit from (−2,0,0) gives
t=1 at (−1,0,0) with normal (−1,0,0); an interior ray exits at +x; parallel and
offset rays miss. WASM smoke checks the radial hit.

## Lean target

Not yet formalized.
