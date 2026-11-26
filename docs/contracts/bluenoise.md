# Contract: blue-noise surface sampling

**Status:** implemented (Wave 2).

## Source (paper only)

- Robert Bridson, "Fast Poisson Disk Sampling in Arbitrary Dimensions",
  SIGGRAPH 2007 Technical Sketch. Adapted to triangle meshes with
  area-weighted triangle selection; the spatial-hash acceleration is a future
  optimization (current implementation uses brute-force separation checks).

## API / constraints

`geomlib_blue_noise_surface(mesh, min_distance, rng, out)` returns positions,
normals and source triangle IDs. A NULL RNG uses a fixed deterministic seed.
C11, wasm32 safe.

## Acceptance

`tests/test_bluenoise.cpp`: 572 samples on a unit sphere at min_distance 0.12;
separation respected; deterministic for a fixed seed. WASM smoke checks the
returned point array.

## Lean target

Not applicable (sampling algorithm).
