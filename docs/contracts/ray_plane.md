# Contract: ray/plane intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard plane-equation intersection: `t = n·(p0−o) / n·d`.

## Approach

Solve for the ray parameter `t` along a direction against an infinite plane.
Reject parallel rays and negative `t`; the segment variant additionally
rejects `t` outside [0,1].

## API / constraints

`geomlib_ray_plane_intersection(origin, direction, plane_point, plane_normal,
t_out, point_out)` and `geomlib_segment_plane_intersection(a, b, ...)`. The
plane normal need not be unit. C11, wasm32 safe.

## Acceptance

`tests/test_ray_plane.cpp`: a ray along +z hits plane z=2 at t=2; away-facing
and parallel rays miss; a segment to z=4 crosses at t=0.5 and a short segment
misses. WASM smoke checks a z=3 hit.

## Lean target

Not yet formalized.
