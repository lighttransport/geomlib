# Contract: ray–cone intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Analytic quadratic intersection of a ray with a right circular cone.

## Approach

Express the cone equation `(p . w)^2 = |p|^2 cos^2(theta)` with `p` relative to
the apex, solve the resulting quadratic in the ray parameter, and keep the
nearest non-negative root lying on the axis-facing nappe.

## API / constraints

`geomlib_ray_cone_intersect(origin, direction, apex, axis, half_angle, out)`
accepts a half-angle in (0, pi/2) and any non-zero direction. C11, wasm32 safe.

## Acceptance

`tests/test_ray_cone.cpp`: a 45-degree cone at the origin hit from (−3,0,2)
along +x gives t=1 at (−2,0,2); an axial ray misses; a parallel ray at x=3 hits
at z=3. WASM smoke checks the side hit.

## Lean target

Not yet formalized.
