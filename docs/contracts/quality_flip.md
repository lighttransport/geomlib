# Contract: surface mesh edge-flip quality improvement

**Status:** implemented (Wave 3).

## Source (paper only)

- Frey and George, "Mesh Generation: Application to Finite Elements", 2000
  (edge-flip improvement of triangle quality).

## Approach

Group triangle edges by endpoint pair; interior edges have two incident faces.
For each, project the quadrilateral onto the plane of one face, require the two
diagonals to cross, flip when the minimum corner angle improves, and orient the
new faces consistently with the old normals. The edge list is rebuilt after
each accepted flip to avoid stale adjacency.

## API / constraints

`geomlib_improve_mesh_quality(mesh, max_passes)` returns the flip count and
modifies triangle indices in place; `geomlib_mesh_min_angle` reports quality.
Boundary edges are never flipped. C11, wasm32 safe.

## Acceptance

`tests/test_quality_flip.cpp`: a distorted 5×5 planar grid improves minimum
angle from 0.0515 to 0.1754 radians in 10 flips while preserving area (16.0000),
boundary edge set and triangle count. WASM smoke checks a positive flip count.

## Lean target

Not yet formalized.
