# Contract: vertex-clustering decimation

**Status:** implemented (Wave 3).

## Source (paper only)

- Low & Tan, "Model Simplification Using Vertex-Clustering", 1997.

## Approach

Hash vertices into uniform voxel cells, average each cluster's positions,
remap triangle indices to cluster representatives, and drop degenerate and
duplicate triangles. Unlike QEM this can change topology and may open holes,
which is acceptable for coarse approximations and collision proxies.

## API / constraints

`geomlib_decimate_vertex_clustering(input, voxel_size, output)`.
C11, wasm32 safe.

## Acceptance

`tests/test_clustering.cpp`: a 1024-triangle sphere simplifies to 747 triangles
with valid in-range indices and finite positions. WASM smoke checks success.

## Lean target

Not yet formalized.
