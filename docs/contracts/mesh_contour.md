# Contract: mesh contour lines

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard marching-triangles applied to per-vertex scalar values on a
  triangle mesh.

## Approach

For each triangle, linearly interpolate edge crossings of the scalar field at
the requested iso value and join the two crossings into a segment. Ambiguous
four-crossing cases are connected by minimum-distance pairing.

## API / constraints

`geomlib_mesh_contour_lines(mesh, values, iso, out)` and
`geomlib_mesh_contour_free`. Values must be finite and one per mesh vertex.
C11, wasm32 safe.

## Acceptance

`tests/test_mesh_contour.cpp`: a 5×5 grid with field `x` at iso 2 gives 8
segments, all endpoints at x=2 with zero value error and vertical orientation;
an out-of-range iso gives no segments. WASM smoke checks a positive count.

## Lean target

Not yet formalized.
