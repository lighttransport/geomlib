# Contract: 2D straight skeleton

**Status:** implemented (Wave 3).

## Sources (papers only)

- Aichholzer, Alberts, Aurenhammer, Krasser, "A Straight Skeleton of a Simple
  Polygon", 1995.
- Felkel & Obdrzalek, "Straight Skeleton Implementation", 1998
  (wavefront propagation approach).

## Approach

Represent the inward-offset wavefront as polygon vertices moving along their
angle bisectors at the speed corresponding to unit-speed edge offsets. At
small time steps, collapse edges that have become degenerate and trace each
surviving vertex's path. The traced paths form the skeleton graph with a
collapse time per node.

## API / constraints

`geomlib_compute_straight_skeleton(polygon, out)` for simple CCW 2D polygons.
C11, wasm32 safe.

## Acceptance

`tests/test_straight_skeleton.cpp`: a unit square traces 202 nodes / 198 edges,
all finite and in bounds, with at least one point near the centre. WASM smoke
passes a square as a flat coordinate array.

## Lean target

Not yet formalized.
