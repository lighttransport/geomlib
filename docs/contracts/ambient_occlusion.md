# Contract: mesh ambient occlusion

**Status:** implemented (Wave 3).

## Source (paper only)

- Zhukov et al., "An Ambient Occlusion Approach to Shading of 3D Scenes",
  GraphiCon 1998 (hemisphere ray sampling).

## Approach

For each vertex, sample a deterministic Fibonacci hemisphere in the tangent
frame of its normal, cast occlusion rays through a BVH, and report the fraction
of unblocked directions. Normals are computed internally, so the mesh argument
is non-const.

## API / constraints

`geomlib_mesh_ambient_occlusion(mesh, samples, max_distance, out)` writes one
value in [0,1] per vertex (1 = fully open). `samples` is clamped to [1, 2048].
C11, wasm32 safe.

## Acceptance

`tests/test_ambient_occlusion.cpp`: a flat grid gives AO 1.0000 everywhere; a
three-wall inside corner gives 0.9531 at the corner vertex. WASM smoke checks
the plane is fully open.

## Lean target

Not yet formalized.
