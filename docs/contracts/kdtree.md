# Contract: k-d tree for point sets

**Status:** implemented (Wave 3).

## Source (paper only)

- Bentley, "Multidimensional Binary Search Trees Used for Associative
  Searching", CACM 1975.

## Approach

Build a balanced tree by recursively splitting the widest coordinate axis at
the median. kNN uses branch-and-bound with a bounded max-heap of size k;
radius search prunes subtrees whose splitting plane is farther than the query
radius.

## API / constraints

`geomlib_kdtree_build`, `geomlib_kdtree_knn`, `geomlib_kdtree_radius`,
`geomlib_kdtree_free`. The tree borrows the point array. C11, wasm32 safe.

## Acceptance

`tests/test_kdtree.cpp`: kNN distances sorted and nearest matches brute force;
radius query count matches brute force (61 points). WASM smoke queries four
points.

## Lean target

Not yet formalized.
