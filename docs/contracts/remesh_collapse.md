# Contract: isotropic edge-collapse coarsening

**Status:** implemented (Wave 3).

## Source (paper only)

- Hoppe, DeRose, Duchamp, McDonald, Stuetzle, "Mesh Optimization", SIGGRAPH
  1993 (edge collapse with link constraints).

## Approach

Find the shortest interior edge above numerical zero, require the link
condition (common neighbours are exactly the two opposite vertices, which must
not already be connected), collapse one endpoint to the midpoint, remove
degenerate and duplicate faces, and compact unreferenced vertices.

## API / constraints

`geomlib_coarsen_edges(mesh, min_edge_length, max_passes, preserve_boundary)`
modifies the mesh in place and returns the collapse count. Each pass performs
one collapse. C11, wasm32 safe.

## Acceptance

`tests/test_remesh_collapse.cpp`: a refined sphere coarsens from 936 to 558
triangles and 497 to 308 vertices in 189 collapses, stays manifold, and keeps
volume within 4%. WASM smoke checks a positive collapse count.

## Lean target

Not yet formalized.
