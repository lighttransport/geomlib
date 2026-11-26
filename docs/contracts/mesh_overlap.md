# Contract: broad-phase AABB overlap pairs

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard bounding-volume-hierarchy broad phase for collision detection; see
  the pbr-book BVH chapter and Ericson, "Real-Time Collision Detection".

## Approach

Build a binary box tree over each mesh by recursively splitting triangle index
ranges in half, accumulating child bounds. Enumerate overlapping pairs by
recursive descent, pruning nodes whose AABBs do not intersect. Leaves fall back
to per-triangle AABB tests.

## API / constraints

`geomlib_aabb_overlap_pairs(mesh_a, mesh_b, out)` fills an owned
`geomlib_overlap_pairs` with two triangle indices per overlapping pair. The
query is broad-phase only: overlapping AABBs do not imply intersecting
triangles. C11, wasm32 safe.

## Acceptance

`tests/test_mesh_overlap.cpp`: identical boxes give 120 overlapping pairs (all
indices in range), boxes 5 apart give 0, boxes shifted by 0.5 give 80. WASM
smoke checks that distant boxes are culled to 0 pairs.

## Lean target

Not yet formalized.
