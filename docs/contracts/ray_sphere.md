# Contract: ray–sphere intersection

**Status:** implemented (Wave 3).

## Source (paper only)

- Analytic quadratic intersection of a ray with a sphere in 3D Euclidean space.

## Approach

Substitute the parametric ray into the sphere equation, yielding
`a t² + b t + c = 0`. Solve via the discriminant; report both roots, entry and
exit points and their outward normals. The direction need not be normalized.

## API / constraints

`geomlib_ray_sphere_intersect(origin, direction, center, radius, out)` returns
1 on intersection and fills `geomlib_ray_sphere_hit`. Negative roots indicate
the sphere lies partially behind the ray origin. C11, wasm32 safe.

## Acceptance

`tests/test_ray_sphere.cpp`: a ray from the origin down +x to a unit sphere at
x=5 gives roots 4 and 6 with correct normals; a perpendicular ray misses; a
tangent ray has coincident roots; an interior origin gives roots −1 and +1.
WASM smoke checks the entry root.

## Lean target

Not yet formalized.
