# Contract: mesh-mesh minimum distance

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard BVH closest-point query; intersecting surfaces are detected with
  the existing mesh-intersection primitive.

## Approach

Build BVHs for both meshes. Query each vertex of A against B and each vertex of
B against A, keep the smallest distance and its point pair. If the meshes
intersect, report zero distance directly.

## API / constraints

`geomlib_mesh_mesh_distance(a, b, out)` returns closest points in world
coordinates. Vertex sampling can slightly overestimate the true surface
distance for nearly touching non-intersecting meshes. C11, wasm32 safe.

## Acceptance

`tests/test_mesh_distance.cpp`: two unit spheres with centres 3 apart give
distance 1.0000 with closest points facing each other; shifted to overlap, the
distance is exactly 0. WASM smoke checks the gap.

## Lean target

Not yet formalized.
