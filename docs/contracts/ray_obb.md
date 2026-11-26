# Contract: ray–OBB intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Slab intersection method in the oriented box's local frame; see Kay and
  Kajiya, "Ray Tracing Complex Scenes", SIGGRAPH 1986.

## Approach

Transform the ray into the OBB's orthonormal local frame, perform the standard
slab test against the half-extents, and return the nearest non-negative root
with the world-space outward face normal. Rays starting inside exit through the
far slab.

## API / constraints

`geomlib_ray_obb_intersect(origin, direction, obb, out)` uses the existing
`geomlib_obb` type. C11, wasm32 safe.

## Acceptance

`tests/test_ray_obb.cpp`: an axis-aligned unit box hit from (−2,0,0) gives t=1
at (−1,0,0) with normal (−1,0,0); a rotated box is hit on the correct local
face; an interior ray exits outward; a ray pointing away misses. WASM smoke
checks the axis-aligned hit.

## Lean target

Not yet formalized.
