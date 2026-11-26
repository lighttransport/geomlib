# Contract: mean-curvature-flow curve skeleton

**Status:** implemented (Wave 3).

## Source (paper only)

- Dey, Sun, Cohen, "Mean Curvature Flow of 3D Meshes via Curvature-driven Edge
  Collapse", IEEE TVCG 2007.

## Approach

Copy the input mesh; repeatedly collapse the shortest edge whose incident face
count is manifold (1 or 2), moving the surviving vertex to the edge midpoint
and compacting unreferenced vertices. Track a representative original vertex
per surviving node. The contracted mesh's vertices and edges form the 1D
skeleton graph. The input mesh is not modified.

## API / constraints

`geomlib_extract_curve_skeleton(mesh, params, out)` with `geomlib_curve_skeleton`
holding points, edges and original vertex representatives. C11, wasm32 safe.

## Acceptance

`tests/test_skeleton.cpp`: a 231-vertex sphere contracts to 131 nodes / 357
edges with valid finite positions and edge indices; input unchanged. WASM smoke
checks node and edge counts.

## Lean target

Topology preservation under the link condition is future work (related to
`formal/GeomLibFormal/Topology.lean`).
