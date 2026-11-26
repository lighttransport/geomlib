# Contract: marching tetrahedra

**Status:** implemented (Wave 3).

## Sources (papers only)

- Gallagher, "Marching Tetrahedra: A New Method for Surface Extraction",
  technical report (1991); standard 6-tetrahedra cube decomposition.

## Approach

Split each grid cube into six tetrahedra around the body diagonal 0->6. For
each tetrahedron, find the crossing edges of the zero level set (3 or 4) and
emit one or two triangles. Triangle winding is fixed by comparing the
inside/outside corner centroids; the 4-crossing quad is ordered into a cycle
so neighbouring tetrahedra agree. Resulting triangles are welded with the
existing soup-to-mesh pipeline.

## API / constraints

`geomlib_sdf_marching_tets(sdf, out)` for dense SDFs. C11, wasm32 safe.

## Acceptance

`tests/test_marching_tets.cpp`: unit sphere at resolution 32 gives 28688
triangles, volume 4.1048 (exact 4.1888) and is watertight. WASM smoke extracts
a sphere surface.

## Lean target

Marching-cubes combinatorics are formalized in `formal/GeomLibFormal/MarchingCubes.lean`;
the tetrahedral table is future work.
