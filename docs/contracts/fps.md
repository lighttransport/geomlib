# Contract: farthest-point sampling

**Status:** implemented (Wave 3).

## Source (paper only)

- Eldar, Lindenbaum, Porat, "The Farthest Point Traversal Algorithm for
  Cost-Efficient Sampling", IEEE TPAMI 2005.
- Moosmann, Nowicki, Schmid, "Random sampling for fast visual object
  recognition", CVPR 2007 (farthest-point sampling on features).

## Approach

Start from a seed vertex; at each step compute the geodesic distance from every
surface vertex to the nearest sample and add the farthest vertex. Radii are
recorded per step, with the final entry equal to the coverage after all samples
have been added. Ties resolve to the lowest vertex index. The heat-method
solver falls back to discrete edge propagation when the continuous solve is
degenerate.

## API / constraints

`geomlib_farthest_point_sampling(mesh, k, seed, out)` and
`geomlib_fps_result_free`. C11, wasm32 safe.

## Acceptance

`tests/test_fps.cpp`: 8 samples on a sphere produce radii
0, 3.941, 3.159, …, 1.267, with distinct in-range vertices, monotone radii,
exact final geodesic coverage, and deterministic repeats. WASM smoke checks 6
samples.

## Lean target

Not yet formalized.
