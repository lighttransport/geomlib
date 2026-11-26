# Contract: variational shape approximation

**Status:** implemented (Wave 3).

## Source (paper only)

- Cohen-Sfady, Alliez, Desbrun, "Variational Shape Approximation", ACM
  SIGGRAPH 2004.

## Approach

Farthest-point seeding of planar proxies, then alternating cluster assignment
(minimizing the sum of squared point-to-plane distances of each triangle's
vertices) and weighted PCA plane refitting. The coarse mesh is a fan from
each proxy point to its region's boundary loops.

## API / constraints

`geomlib_vsa_simplify(input, region_count, iterations, result, output)` returns
both the triangle clustering with proxy planes and a simplified triangle mesh.
C11, wasm32 safe.

## Acceptance

`tests/test_vsa.cpp`: a 576-triangle sphere approximated by 12 regions produces
a 38-triangle mesh with valid region ids. WASM smoke checks an 8-region run.

## Lean target

Not yet formalized.
