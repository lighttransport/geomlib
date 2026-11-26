# Contract: 2D polygon rasterization

**Status:** implemented (Wave 3).

## Source (paper only)

- Cell-centre point-in-polygon sampling with signed boundary distance.

## Approach

Compute the polygon bounds with a margin, create a uniform grid covering them,
and test each cell centre for containment. Optional signed distances are the
minimum distance to the polygon edges, negative inside.

## API / constraints

`geomlib_rasterize_polygon2d(polygon, n, nx, ny, margin, with_distance, out)`
fills an owned `geomlib_raster2d`. C11, wasm32 safe.

## Acceptance

`tests/test_raster2d.cpp`: a 2×2 square at 50×50 gives area 4.0643 (analytic 4)
and signed distances consistent with the mask. WASM smoke checks the inside count.

## Lean target

Not yet formalized.
