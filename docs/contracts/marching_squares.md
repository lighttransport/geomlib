# Contract: marching squares

**Status:** implemented (Wave 3).

## Source (paper only)

-lorensen & Cline, "Marching Cubes", SIGGRAPH 1987 (2D analogue).

## Approach

Classify each grid cell's four corners, emit one or two line segments from the
zero crossings on cell edges, use the cell-centre average as an asymptotic
decider for saddle cells, then chain segments into polylines by welding shared
endpoints.

## API / constraints

`geomlib_marching_squares(values, nx, ny, h, origin, out)` for node-sampled
scalar grids. C11, wasm32 safe.

## Acceptance

`tests/test_marching_squares.cpp`: a 64x64 circle field yields one contour with
perimeter 6.3319 vs analytic 6.2832. WASM smoke extracts a circle contour.

## Lean target

Case-table combinatorics are partially formalized in
`formal/GeomLibFormal/MarchingCubes.lean`; the 2D table is future work.
