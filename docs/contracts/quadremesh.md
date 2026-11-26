# Contract: quad remeshing

**Status:** implemented (Wave 3).

## Source (paper only)

- Botsch, Kobbelt, Pauly, Alliez, Levy, "Polygon Mesh Processing", AK Peters
  2010, chapter on remeshing (triangle-to-quad conversion by edge pairing).

## Approach

Build triangle adjacency from sorted undirected edges; greedily merge adjacent
triangle pairs into quads when the dihedral angle and all four interior quad
angles satisfy bounds. Unmerged triangles are kept. An optional relaxation
moves shared vertices toward the average of their incident quad centroids.

## API / constraints

`geomlib_remesh_to_quads(input, params, geomlib_quadmesh*)`; the output is a
mixed quad/triangle mesh with shared vertices and explicit `quadmesh_free`.
C11, wasm32 safe.

## Acceptance

`tests/test_quadremesh.cpp`: a planar 5x5 grid yields 21 quads and 8 triangles,
element count conserved, all quads have four distinct in-range vertices.
WASM smoke checks a nonzero quad count.

## Lean target

Quad validity angle bounds are future work.
