# Contract: triangle-mesh intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Moller and Trumbore, "Fast, Minimum Storage Ray/Triangle Intersection", 1997.
- Standard face-pair classification of intersection polylines.

## Approach

Trace every edge of A through a BVH of B and every edge of B through a BVH of
A. Each crossing is tagged with its (A face, B face) pair; points are grouped by
that pair, ordered by angle around their centroid in the A face plane, and
joined into segments. Crossings must lie within both AABBs, which rejects
spurious tangent hits.

## API / constraints

`geomlib_mesh_mesh_intersection(a, b)` and `geomlib_intersection_free`.
Coplanar overlapping faces are not reported; the angular ordering assumes a
convex cross-section per face pair. C11, wasm32 safe.

## Acceptance

`tests/test_intersection.cpp`: two overlapping unit boxes (offset 0.5/0.3/0.2)
give 19 segments, all endpoints inside the overlap and within 1e-3 of both
surfaces; disjoint boxes give 0. WASM smoke checks a non-empty result.

## Lean target

Not yet formalized.
