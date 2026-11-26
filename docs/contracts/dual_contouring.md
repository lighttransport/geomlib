# Contract: dual contouring

**Status:** implemented (Wave 3).

## Source (paper only)

- Ju, Losasso, Schaefer, Hanrahan, "Dual Contouring of Hermite Data",
  SIGGRAPH 2002.

## Approach

For each crossing grid cell, gather Hermite data (edge intersection points and
SDF gradients), build a quadric per point and solve a 3x3 system for a single
representative vertex (falling back to the average of intersections when
singular). Around each crossing grid edge, emit two triangles connecting the
four surrounding cell vertices. Winding is fixed with the SDF gradient. The
resulting soup is welded with the existing soup-to-mesh pipeline.

## API / constraints

`geomlib_sdf_dual_contouring(sdf, out)` for dense SDFs. C11, wasm32 safe.

## Acceptance

`tests/test_dual_contouring.cpp`: unit sphere at resolution 24 gives 5376
triangles, volume 4.0886 (exact 4.1888) and is watertight. WASM smoke checks
a smaller sphere.

## Lean target

Not yet formalized.
