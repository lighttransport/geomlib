# Contract: 2D polygon triangulation with holes

**Status:** implemented (Wave 3).

## Source (paper only)

- Chew, "Constrained Delaunay Triangulations", Algorithmica 1989 (CDT is used
  to respect all boundary edges).
- Standard hole removal by centroid classification.

## Approach

Concatenate outer and hole vertices, constrain every boundary edge with CDT,
then discard triangles whose centroid lies inside a hole. Triangle indices
refer to the concatenated vertex array.

## API / constraints

`geomlib_triangulate_with_holes(outer, no, holes, hole_sizes, nh, out)` and
`geomlib_polytri_result_free`. Holes must lie strictly inside the outer polygon
and not overlap each other. C11, wasm32 safe.

## Acceptance

`tests/test_polytri_holes.cpp`: a 4×4 square with a 2×2 hole yields 8 triangles
of total area 12.0000 with none covering the hole; an 8×4 box with two 1×1
holes yields area 30.0000. WASM smoke checks area 3 for a 2×2 square with a 1×1
hole.

## Lean target

Not yet formalized.
