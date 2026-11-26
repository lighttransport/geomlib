# Contract: 2D polygon offsetting

**Status:** implemented (Wave 3).

## Sources (papers only)

- Preparata & Shamos, "Computational Geometry: An Introduction", 1985
  (line-offset construction); standard miter join treatment.

## Approach

Translate each polygon edge along its inward normal by the signed distance
(positive distance offsets outward), intersect consecutive offset lines to
obtain new vertices, clamp over-long miters, and split self-crossing offset
loops at their intersections into simple output loops.

## API / constraints

`geomlib_offset_polygon(polygon, distance, arc_segments, out)`; `arc_segments`
is reserved for future round joins. C11, wasm32 safe.

## Acceptance

`tests/test_polygon_offset.cpp`: unit square offsets to area 1.44 (outward
0.1) and 0.64 (inward 0.1), exactly. WASM smoke checks a single 4-vertex loop.

## Lean target

Not yet formalized.
