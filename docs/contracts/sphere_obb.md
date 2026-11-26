# Contract: sphere–OBB overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Local-frame closest-point test between a sphere and an oriented box.

## Approach

Express the sphere centre in the OBB's orthonormal frame, clamp each coordinate
to the half-extents, reconstruct the world-space closest point, and compare the
centre-to-box distance with the radius.

## API / constraints

`geomlib_sphere_obb_overlap(center, radius, obb, out)` uses the existing
`geomlib_obb` type. C11, wasm32 safe.

## Acceptance

`tests/test_sphere_obb.cpp`: a unit half-extent box hit from (2,0,0) misses at
radius 0.5 and penetrates by 0.5 at radius 1.5; a rotated box gives the same
gap in local coordinates. WASM smoke checks the penetrating case.

## Lean target

Not yet formalized.
