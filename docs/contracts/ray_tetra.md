# Contract: ray–tetrahedron intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Möller-Trumbore face intersection applied to the four faces of a tetrahedron.

## Approach

Test the ray against all four faces with the shared ray/triangle primitive,
select the nearest non-negative hit and convert the normalized-ray parameter
back to the caller's direction units. The face normal follows the face winding.

## API / constraints

`geomlib_ray_tetrahedron_intersect(origin, direction, v0, v1, v2, v3, out)`
accepts any non-zero direction. C11, wasm32 safe.

## Acceptance

`tests/test_ray_tetra.cpp`: a ray from (−1,−1,−1) along (1,1,1) enters at the
origin (t=1) inside the solid; a ray from inside exits; a ray pointing away
misses. WASM smoke checks the entry parameter.

## Lean target

Not yet formalized.
