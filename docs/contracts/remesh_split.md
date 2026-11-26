# Contract: isotropic edge-split refinement

**Status:** implemented (Wave 3).

## Source (paper only)

- Frey and George, "Mesh Generation: Application to Finite Elements", 2000
  (edge-split refinement for triangular surfaces).

## Approach

Find the longest edge exceeding the target, insert its midpoint as a new
vertex, and replace its one incident triangle with two (boundary edge) or two
incident triangles with four (interior edge). New faces are oriented against
the original face normal. Adjacency is rebuilt from scratch each pass.

## API / constraints

`geomlib_refine_edges(mesh, max_edge_length, max_passes)` modifies the mesh in
place and returns the split count. Each pass performs one split, so `max_passes`
bounds the number of inserted vertices. Only refinement is performed; there is
no edge collapse for coarsening. C11, wasm32 safe.

## Acceptance

`tests/test_remesh_split.cpp`: a coarse unit sphere refined to 0.35 gives 302
splits, longest edge 0.8805 → 0.3314, volume 3.5186 → 3.3565, and remains
manifold. WASM smoke checks splits occur and the longest edge shrinks.

## Lean target

Not yet formalized.
