# Contract: bounded Hausdorff distance

**Status:** implemented (Wave 0). **Source:** 1-Lipschitz property of the
unsigned distance function to a closed set (standard analysis fact);
BVH-accelerated closest points from the existing geomlib BVH module.

## Specification

For source mesh S and target T, the unsigned distance d(x,T) is 1-Lipschitz,
so for any triangle with vertices v_i:

    sup_{x in triangle} d(x,T) <= min_i d(v_i,T) + longest_edge(triangle)

Sampling vertices, edge midpoints and centroids gives a lower bound; the
inequality gives a rigorous upper bound. Both directions are combined.

## API / constraints

C11, uses BVH; no deep recursion. WASM: yes.

## Acceptance

`tests/test_distance.cpp`: lower ≤ upper; identical meshes have zero sampled
distance; a 0.25 translation produces lower bound 0.25 and an upper bound that
contains the truth.

## Lean target

The Lipschitz inequality can be formalized in a future analysis layer; not yet
proved.
