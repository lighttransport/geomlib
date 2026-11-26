# Contract: 2D Delaunay refinement

**Status:** implemented (Wave 3).

## Source (paper only)

- Chew, "Guaranteed-Quality Triangular Meshes", 1989.
- Ruppert, "A Delaunay Refinement Algorithm for Quality 2-Dimensional Mesh
  Generation", 1995.

## Approach

Repeatedly build a constrained Delaunay triangulation of the current points and
boundary constraints. Insert the circumcentre of the worst triangle when it is
strictly inside the domain; otherwise split that triangle's longest edge at its
midpoint and update the constraint list. Rebuild until the angle target is met
or the Steiner-point budget runs out.

## API / constraints

`geomlib_refine_polygon2d(polygon, n, params, out)` with
`geomlib_refine2d_params_default` and `geomlib_cdt2d_min_angle`. The result owns
its vertices and triangles. This is a best-effort refinement: highly elongated
domains can retain self-similar skinny triangles, so the angle target is not a
formal guarantee. C11, wasm32 safe.

## Acceptance

`tests/test_refine2d.cpp`: a 4×1 rectangle gains Steiner points (24 vertices,
34 triangles) while area stays 4.0000, all triangles remain inside the domain,
and the minimum angle does not get worse. WASM smoke checks Steiner insertion.

## Lean target

Not yet formalized.
