# Contract: circle intersection and union area

**Status:** implemented (Wave 3).

## Source (paper only)

- Analytic lens area formula for two intersecting circles.

## Approach

Let the radii be a, b and centre distance d. For disjoint or tangent circles
return 0; for containment return `pi*min(a,b)^2`; otherwise use
`a^2 acos(c0) + b^2 acos(c1) - 0.5*sqrt(Heron(a,b,d))` with the cosine-law
arguments. Union area is the sum of disc areas minus the lens.

## API / constraints

`geomlib_circle_intersection_area(r0,r1,d)` and
`geomlib_circle_union_area(r0,r1,d)` accept absolute radii and distance. C11,
wasm32 safe.

## Acceptance

`tests/test_circle_area.cpp`: coincident unit circles give pi; tangent and
disjoint circles give 0; d=1 gives `2pi/3 - sqrt(3)/2 = 1.228370`; the union
matches the sum minus the lens. WASM smoke checks the lens value.

## Lean target

Not yet formalized.
