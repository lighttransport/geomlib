# Contract: plane-plane intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Analytic intersection of two planes in 3D Euclidean space.

## Approach

Normalize both plane normals; their cross product is the line direction. Solve
the two plane equations for a point on the line with the cross-product formula
`p = ((d0*n1 - d1*n0) x dir) / |dir|^2`. Parallel planes return a `parallel`
flag.

## API / constraints

`geomlib_plane_plane_intersection(point0, normal0, point1, normal1, out)` uses
the same plane convention as `geomlib_ray_plane_intersection`. C11, wasm32 safe.

## Acceptance

`tests/test_plane_plane.cpp`: z=0 and x=1 intersect in a unit-direction y-axis
line through (1,0,0); points along it satisfy both planes; parallel planes
report no line. WASM smoke checks the same perpendicular case.

## Lean target

Not yet formalized.
