# Contract: triangle–AABB overlap

**Status:** implemented (Wave 3).

## Source (paper only)

- Akenine-Möller, "Fast 3D Triangle-Box Overlap Testing", 2001; the 13-axis
  separating-axis test.

## Approach

Translate the triangle so the box is centred at the origin, then test separation
along the three box face normals, the triangle normal and the nine cross
products of box axes with triangle edges. A generic projection routine handles
all axes uniformly.

## API / constraints

`geomlib_triangle_aabb_overlap(v0,v1,v2,bmin,bmax)` returns 1 on overlap.
Degenerate axes are skipped. C11, wasm32 safe.

## Acceptance

`tests/test_tri_aabb.cpp`: a triangle crossing the unit box and a triangle
inside it overlap; triangles above or beside the box are separated. WASM smoke
checks the crossing case.

## Lean target

Not yet formalized.
