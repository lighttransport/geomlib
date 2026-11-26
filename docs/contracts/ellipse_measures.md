# Contract: ellipse area and perimeter

**Status:** implemented (Wave 3).

## Source (paper only)

- Ramanujan, "Modular Equations and Approximations to pi", Quarterly Journal of
  Mathematics 1914 (second perimeter approximation).

## Approach

Area is `pi*a*b`. The perimeter uses
`C ≈ pi*(a+b)*(1 + 3h/(10+sqrt(4-3h)))` with `h = ((a-b)/(a+b))^2`, which is
exact for circles and highly accurate for moderate eccentricities.

## API / constraints

`geomlib_ellipse_area(a,b)` and `geomlib_ellipse_perimeter(a,b)` accept absolute
semi-axis lengths. C11, wasm32 safe.

## Acceptance

`tests/test_ellipse_measures.cpp`: a unit circle gives 2pi exactly; semi-axes
(2,1) give perimeter 9.688448 (true ≈ 9.68845) and area 2pi. WASM smoke checks
the perimeter.

## Lean target

Not yet formalized.
