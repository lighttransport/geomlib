# Contract: 2D polygon Poisson-disk sampling

**Status:** implemented (Wave 3).

## Source (paper only)

- Bridson, "Fast Poisson Disk Sampling in Arbitrary Dimensions", SIGGRAPH
  Sketches 2007.

## Approach

Bridson dart throwing in a background grid, with each candidate accepted only
if it lies inside the simple polygon (ray-crossing test) and is at least `radius`
from neighbouring grid cells. Deterministic PCG32 seeding.

## API / constraints

`geomlib_poisson_disk_polygon(polygon, np, params, out, capacity)` returns the
sample count. C11, wasm32 safe.

## Acceptance

`tests/test_poisson2d.cpp`: a 2×2 box with r=0.25 gives 46 samples, minimum
separation 0.2511, all inside and deterministic. WASM smoke checks a unit box.

## Lean target

Not yet formalized.
