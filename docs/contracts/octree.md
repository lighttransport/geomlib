# Contract: point-cloud octree

**Status:** implemented (Wave 3).

## Source (paper only)

- Meagher, "Geometric Modeling Using Octree Encoding", Computer Graphics and
  Image Processing, 1982.

## Approach

Build an adaptive octree by recursively splitting overloaded nodes into eight
octants. Radius queries prune nodes whose AABB is farther than the query
radius; kNN uses best-first node traversal with an ordered candidate buffer.

## API / constraints

`geomlib_octree_build`, `geomlib_octree_radius`, `geomlib_octree_knn`,
`geomlib_octree_free`. The tree borrows the point array. C11, wasm32 safe.

## Acceptance

`tests/test_octree.cpp`: 240-point cloud gives nearest neighbour 0.0874
matching brute force and radius count 112 matching brute force. WASM smoke
runs a radius query.

## Lean target

Not yet formalized.
