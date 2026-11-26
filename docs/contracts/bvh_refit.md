# Contract: dynamic BVH refit

**Status:** implemented (Wave 3).

## Source (paper only)

- Typical BVH refit strategy (Wald et al.): recompute leaf bounds and union
  child bounds bottom-up without rebuilding topology.

## Approach

For each leaf, recompute its AABB from the current positions of its triangles
(using the BVH's existing primitive order). Interior bounds are unions of
child bounds, computed by post-order traversal. The mesh pointer is updated.

## API / constraints

`geomlib_bvh_refit(bvh, mesh)`; requires identical mesh topology to the build.
C11, wasm32 safe.

## Acceptance

`tests/test_bvh_refit.cpp`: after translating a sphere by (3,0,0), refitted
closest-point distance matches a freshly rebuilt tree (0.5 vs 0.5).
WASM smoke checks a refitted query.

## Lean target

Not yet formalized.
