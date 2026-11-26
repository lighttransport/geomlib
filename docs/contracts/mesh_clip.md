# Contract: mesh plane clipping

**Status:** implemented (Wave 3).

## Source (paper only)

- Standard half-space clipping of triangles with cross-section cap reconstruction.

## Approach

Classify each vertex against the plane, clip each triangle to the kept
half-space (producing 3 or 4 vertices), fan-triangulate the clipped polygon,
then collect cut-edge points, deduplicate, sort by angle in the plane and fan-
triangulate the cap with the orientation chosen to close the solid.

## API / constraints

`geomlib_clip_mesh_plane(mesh, point, normal, keep_positive, out)` assumes a
single connected cross-section loop (convex cuts). C11, wasm32 safe.

## Acceptance

`tests/test_mesh_clip.cpp`: clipping a unit box by x=0 keeping x<0 gives volume
0.5000. WASM smoke checks the same volume.

## Lean target

Not yet formalized.
