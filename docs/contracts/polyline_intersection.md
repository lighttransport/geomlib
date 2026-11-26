# Contract: polyline / triangle-mesh intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Möller and Trumbore, "Fast, Minimum Storage Ray/Triangle Intersection", 1997.

## Approach

Build a BVH of the mesh, then trace each polyline segment as a bounded ray,
advancing the origin past each hit until the segment is exhausted. Repeated
points are deduplicated within a small tolerance.

## API / constraints

`geomlib_polyline_mesh_intersection(polyline, n, mesh, out)` and
`geomlib_polyline_hits_free`. C11, wasm32 safe.

## Acceptance

`tests/test_polyline_intersection.cpp`: a line through a unit sphere's centre
gives 2 crossings on the surface; a missing line gives 0; a three-point
polyline through the centre also gives 2. WASM smoke checks 2 hits.

## Lean target

Not yet formalized.
