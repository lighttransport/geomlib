# Contract: Ritter bounding sphere

**Status:** implemented (Wave 3).

## Source (paper only)

- Ritter, "An Efficient Bounding Sphere", Graphics Gems, 1990.

## Approach

Pick the farthest pair among the extreme points on the three axes, seed the
sphere with their midpoint and half-distance, then grow the sphere once per
remaining point that lies outside it.

## API / constraints

`geomlib_bounding_sphere_points(points, n, out)` and
`geomlib_bounding_sphere_mesh(mesh, out)`. The result is a valid enclosure but
not the minimum enclosing sphere. C11, wasm32 safe.

## Acceptance

`tests/test_bounding_sphere.cpp`: 300 points on a unit sphere give radius
1.0372 with all points enclosed; a unit box's eight corners give 1.7321 (the
exact half-diagonal). WASM smoke checks the box radius.

## Lean target

Not yet formalized.
