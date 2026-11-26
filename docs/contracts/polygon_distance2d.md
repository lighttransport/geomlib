# Contract: 2D polygon–polygon distance

**Status:** implemented (Wave 3).

## Source (paper only)

- Closest-feature enumeration for simple polygons.

## Approach

Return 0 when a vertex lies inside the other polygon or edges intersect.
Otherwise enumerate vertex-to-edge projections in both directions, using the
existing segment closest-point solver, and keep the minimum with witnesses.

## API / constraints

`geomlib_polygon_distance2d(a, na, b, nb, out)` requires simple polygons with at
least three vertices. C11, wasm32 safe.

## Acceptance

`tests/test_polygon_distance2d.cpp`: unit squares one apart give distance 1 at
(1,0) and (2,0); overlapping squares give 0; diagonal squares give √2. WASM
smoke checks the unit gap.

## Lean target

Not yet formalized.
