# Contract: dihedral-angle region growing

**Status:** implemented (Wave 3).

## Source (paper only)

- Besl and Jain, "Segmentation through Variable-Order Surface Fitting",
  IEEE TPAMI 1988 (region growing by surface continuity).

## Approach

Build face adjacency from shared edges, then grow connected components by BFS,
accepting neighbouring faces whose normal dot product exceeds the cosine of the
dihedral threshold. Each face receives a region id.

## API / constraints

`geomlib_mesh_region_growing(mesh, dihedral_threshold_deg, labels, region_count)`
with `labels` preallocated to `mesh->triangle_count`. Degenerate faces may
split otherwise smooth regions (e.g. sphere pole fans). C11, wasm32 safe.

## Acceptance

`tests/test_region_growing.cpp`: a unit box at 30° gives 6 face regions; a flat
5×5 grid gives 1 region; all labels are valid. WASM smoke checks 6 box regions.

## Lean target

Not yet formalized.
