# Contract: 3D volume Poisson-disk sampling

**Status:** implemented (Wave 3).

## Source (paper only)

- Bridson, "Fast Poisson Disk Sampling in Arbitrary Dimensions", SIGGRAPH
  Sketches 2007.

## Approach

Maintain an active list and a background grid with cell size r/√3. For each
active sample, try up to `candidates` random points in a spherical shell of
radius [r, 2r], accept those inside the box and at least r from neighbouring
grid cells, then continue from accepted samples. PCG32 with a fixed seed makes
the result deterministic.

## API / constraints

`geomlib_poisson_disk_volume(params, out, capacity)` returns the number of
samples written. C11, wasm32 safe.

## Acceptance

`tests/test_poisson3d.cpp`: a 2³ box with r=0.25 gives 153 samples, minimum
separation 0.2504, all inside the box, and identical output for the same seed.
WASM smoke checks a unit box.

## Lean target

Not yet formalized.
