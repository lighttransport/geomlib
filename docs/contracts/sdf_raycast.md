# Contract: SDF ray casting

**Status:** implemented (Wave 3).

## Source (paper only)

- Hart, "Sphere Tracing: A Geometric Method for the Antialiased Ray Tracing of
  Implicit Surfaces", 1996.

## Approach

Sphere trace along the ray using trilinearly interpolated signed distances,
advancing by `max(|d|, eps)` until `|d| <= eps`. When reached, compute the
surface normal by central differences of the SDF gradient. Rays leaving the
grid report a miss.

## API / constraints

`geomlib_sdf_raycast(sdf, ray, max_distance, max_steps, out)`. Requires a dense
SDF with finite values; sparse/narrow-band grids return FLT_MAX outside the
band. The hit tolerance is one quarter of the voxel size. C11, wasm32 safe.

## Acceptance

`tests/test_sdf_raycast.cpp`: a ray from (1.9,0,0) toward a radius-1 sphere hits
at t=0.9005 with radius 0.9995 and normal +x; a ray from inside exits at
t=0.9946; an off-surface ray misses. WASM smoke checks a hit.

## Lean target

Not yet formalized.
