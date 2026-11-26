# Contract: alpha wrapping

**Status:** implemented (Wave 3).

## Source (paper only)

- CGAL "Alpha wrapping" package documentation / concept description
  (watertight wrapping of arbitrary geometry with an alpha offset).

## Approach

Build an unsigned SDF of the input (open surfaces are supported), dilate it by
the alpha distance, optionally smooth the field, and extract the isosurface.
The grid is padded so the offset slab is closed inside the domain.

## API / constraints

`geomlib_alpha_wrap(input, params, out)`; default alpha is 1% of the bounding
box diagonal, resolution 48. C11, wasm32 safe.

## Acceptance

`tests/test_alpha_wrap.cpp`: an open plane wraps into a 5840-triangle watertight
manifold mesh enclosing a nonzero volume. WASM smoke wraps a plane.

## Lean target

Not yet formalized.
