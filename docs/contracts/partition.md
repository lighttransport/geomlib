# Contract: balanced mesh partition

**Status:** implemented (Wave 3).

## Source (paper only)

- Geogram `mesh_partition` concept description; farthest-point Voronoi seeding
  with deterministic rebalancing.

## Approach

Farthest-point seeds in Euclidean space, nearest-seed Voronoi assignment, then
deterministic rebalancing that moves the vertex closest to an underfull seed
from the most overfull part until sizes are balanced.

## API / constraints

`geomlib_partition_mesh(mesh, part_count, out)` returning per-vertex labels and
per-part sizes. C11, wasm32 safe.

## Acceptance

`tests/test_partition.cpp`: a 231-vertex sphere splits into 4 parts of sizes
58/58/58/57; every vertex labelled. WASM smoke checks the part-size array.

## Lean target

Not yet formalized.
