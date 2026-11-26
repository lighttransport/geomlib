# Contract: mesh slicing

**Status:** implemented (Wave 3).

## Source (paper only)

- CGAL `Polygon_mesh_slicer` package description (plane/mesh intersection
  polylines).

## Approach

Clip every triangle by the plane; triangles with sign changes contribute one
segment between the two edge crossings. Segments are chained by matching
endpoints within a welding tolerance into closed/open polylines.

## API / constraints

`geomlib_slice_mesh(mesh, point, normal, weld_eps, out)` returning concatenated
polyline points with per-polyline start indices and counts. C11, wasm32 safe.

## Acceptance

`tests/test_slicer.cpp`: a unit box sliced by z=0 yields one polyline lying on
the plane with at least four points. WASM smoke checks polyline counts.

## Lean target

Not yet formalized.
