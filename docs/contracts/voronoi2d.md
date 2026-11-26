# Contract: 2D Voronoi diagram

**Status:** implemented (Wave 3).

## Source (paper only)

- Voronoi, "Nouvelles applications des paramètres continus à la théorie des
  formes quadratiques", 1908. Half-plane intersection of perpendicular
  bisectors is the classical primal construction.

## Approach

For each site, start with the caller's bounding box polygon (CCW) and clip it
by the bisector half-plane of the site and every other site using
Sutherland–Hodgman. Bisector of p_i and p_j:
`2(p_j - p_i)·x <= |p_j|² - |p_i|²`. Duplicate sites are ignored. Cells are
convex, CCW, and stored as offsets into a flat vertex array.

## API / constraints

`geomlib_compute_voronoi2d(points, n, bmin, bmax, out)`,
`geomlib_voronoi2d_locate`, `geomlib_voronoi2d_cell_area`,
`geomlib_voronoi2d_create/free`. Clipping is O(n²) per diagram; suitable for
moderate point counts. C11, wasm32 safe.

## Acceptance

`tests/test_voronoi2d.cpp`: 25 sites on a 5×5 grid tile the unit box with total
area 1.000000, all 100 cell vertices are nearest to their own site, and 400
nearest-site queries agree with convex-cell containment. WASM smoke checks a
16-site diagram tiles the box.

## Lean target

Not yet formalized.
