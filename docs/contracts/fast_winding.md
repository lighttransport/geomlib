# Contract: fast winding number

**Status:** implemented (Wave 3).

## Source (paper only)

- Barill, Dickson, Wang, Macklin, Wetzstein, Solomon, "Fast Winding Numbers
  for Soups and Clouds", ACM SIGGRAPH 2018.

## Approach

Build a binary tree of triangle clusters split along their longest centroid
axis. Each node stores its area-weighted centroid and normal moment. A far
query uses the dipole solid-angle approximation; near nodes recurse, and leaves
use the exact solid-angle sum. The tree borrows the source mesh (which must
outlive it).

## API / constraints

`geomlib_winding_tree_build(mesh, leaf_size)`,
`geomlib_fast_winding_number(tree, point)`, `geomlib_winding_tree_free`.
C11, wasm32 safe.

## Acceptance

`tests/test_fast_winding.cpp`: sphere inside/outside winding matches the exact
generalized winding number within 0.05. WASM smoke checks the inside value.

## Lean target

Not yet formalized.
