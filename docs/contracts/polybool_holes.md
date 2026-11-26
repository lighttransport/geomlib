# Contract: polygon Boolean with holes

**Status:** implemented (Wave 3).

## Source (paper only)

- Vatti, "A Generic Solution to Polygon Clipping", 1992; hole bridging is a
  standard technique for reducing multiply-connected polygons to simple ones.

## Approach

Represent each polygon as an outer ring plus hole rings. Bridge each hole into
the outer ring by a visible diagonal traversed in reverse winding, forming a
single self-touching simple polygon. Delegate to `geomlib_polygon_boolean`.
Ambiguous four-crossing configurations around bridge slits are resolved
heuristically, so results are best-effort near the bridges.

## API / constraints

`geomlib_polygon_with_holes_set` builds an owned representation;
`geomlib_polygon_boolean_holes(a, b, op, out)` performs union/intersection/
difference. C11, wasm32 safe.

## Acceptance

`tests/test_polybool_holes.cpp`: a framed 4×4 square (2×2 hole) unioned with a
1×1 square inside the hole gives area 13.0000 exactly (frame 12 plus square 1)
and intersection is empty. WASM smoke checks Boolean completion.

## Lean target

Not yet formalized.
