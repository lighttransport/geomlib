# Contract: hex-dominant meshing

**Status:** implemented (Wave 3).

## Source (paper only)

- MeshCrunch / HEXCRUNCH concepts (Lyon et al.); 5-tet split of a hexahedron.

## Approach

Build edge-to-tetrahedra incidence, find edges shared by exactly five
tetrahedra whose union has eight distinct vertices, order the six ring vertices
around the shared diagonal by angle, and emit one hex cell while consuming all
five tetrahedra. Unmerged tetrahedra are passed through.

## API / constraints

`geomlib_hex_dominant(input, out)` returning `geomlib_hexmesh` with borrowed
vertex positions and owned index arrays. C11, wasm32 safe.

## Acceptance

`tests/test_hexdominant.cpp`: a five-tet fan around one edge merges into one
hex with eight distinct vertices and zero remaining tets. WASM smoke runs a
demo merge.

## Lean target

Not yet formalized.
