# Clean-room algorithm contracts

Every new algorithm in geomlib is recorded here **before** implementation.

## Rule

Implementations are written from published papers and mathematical
specifications only. Source code of CGAL, libigl and Geogram is never read or
copied. Their APIs may be used as *feature descriptions*, never as code
templates. Constants and lookup tables are derived from formulas or generated
from the Lean formalization (`formal/`) rather than transcribed.

## Contract template

1. **Bibliographic source** — paper title, authors, year, venue, URL.
2. **Specification** — inputs, outputs, invariants, degeneracy behaviour.
3. **Derivation notes** — how tables/constants were obtained.
4. **API** — geomlib-style signatures (C11, `extern "C"`, init/free).
5. **Acceptance** — analytic fixtures, invariants, golden vectors.
6. **Lean target** — which formalization covers the core property.
7. **WASM capability** — yes/no and any restrictions.

## Contracts

- [Robust predicates](predicates.md) — done (Wave 0)
- [Sparse linear algebra](sparse.md) — done (Wave 0)
- [Mesh I/O](io.md) — done (Wave 0)
- [Bounded Hausdorff distance](distance.md) — done (Wave 0)
- [Parameterization completion](parameterization.md) — done (Wave 1)
- [Wave 1 gap completion](wave1.md) — atlas, Catmull–Clark, repair, vcache — done
- [Mesh booleans](boolean.md) — done (Wave 2)
- [Rigid ICP registration](registration.md) — done (Wave 2)
- [Blue-noise surface sampling](bluenoise.md) — done (Wave 2)
- [Poisson surface reconstruction](poisson.md) — done (Wave 2)
- [ACVD-style remeshing](acvd.md) — done (Wave 2, isotropic CVT).
- [ARAP surface deformation](deform.md) — done (Wave 3).
- [Manifold harmonics](harmonics.md) — done (Wave 3).
- [RANSAC shape detection](shape_detection.md) — done (Wave 3).
- [Quad remeshing](quadremesh.md) — done (Wave 3).
- [Anisotropic remeshing](anisotropic.md) — done (Wave 3).
- [Marching tetrahedra](marching_tets.md) — done (Wave 3).
- [Dual contouring](dual_contouring.md) — done (Wave 3).
- [Mesh saliency](saliency.md) — done (Wave 3).
- [Curve skeleton](skeleton.md) — done (Wave 3).
- [2D straight skeleton](straight_skeleton.md) — done (Wave 3).
- [2D polygon offsetting](polygon_offset.md) — done (Wave 3).
- [Variational shape approximation](vsa.md) — done (Wave 3).
- [Feature-preserving point-set processing](pointset.md) — done (Wave 3).
- [Alpha wrapping](alpha_wrap.md) — done (Wave 3).
- [Mesh slicing](slicer.md) — done (Wave 3).
- [Approximate convex decomposition](convex_decomposition.md) — done (Wave 3).
- [Implicit mean-curvature fairing](fairing.md) — done (Wave 3).
- [Fast winding number](fast_winding.md) — done (Wave 3).
- [Balanced mesh partition](partition.md) — done (Wave 3).
- [Vertex-clustering decimation](clustering.md) — done (Wave 3).
- [Marching squares](marching_squares.md) — done (Wave 3).
- [2D Delaunay triangulation](delaunay2d.md) — done (Wave 3).
- [Fast-marching geodesics](geodesic_fmm.md) — done (Wave 3).
- [Continuous geodesic unfolding](geodesic_exact.md) — done (Wave 3).
- [Linear blend skinning](skinning.md) — done (Wave 3).
- [Dual-quaternion skinning](dualquat_skinning.md) — done (Wave 3).
- [Cross fields](cross_field.md) — done (Wave 3).
- [Shape Diameter Function](shape_diameter.md) — done (Wave 3).
- [Point-set k-d tree](kdtree.md) — done (Wave 3).
- [Dynamic BVH refit](bvh_refit.md) — done (Wave 3).
- [Point-cloud octree](octree.md) — done (Wave 3).
- [Sparse direct solver](sparse_direct.md) — done (Wave 3).
- [Hex-dominant meshing](hexdominant.md) — done (Wave 3).
- [Variational RBF reconstruction](rbf_reconstruct.md) — done (Wave 3).
- [Grid-based fluid solver](fluid.md) — done (Wave 3).
- [2D Voronoi diagram](voronoi2d.md) — done (Wave 3).
- [Triangle-mesh voxelization](voxelize.md) — done (Wave 3).
- [Spatial sorting (Morton/Hilbert)](spatial_sort.md) — done (Wave 3).
- [Triangle-mesh intersection](intersection.md) — done (Wave 3).
- [2D constrained Delaunay](cdt2d.md) — done (Wave 3).
- [SDF-based mesh offset](offset.md) — done (Wave 3).
- [Farthest-point sampling](fps.md) — done (Wave 3).
- [2D polygon medial axis](medial2d.md) — done (Wave 3).
- [2D polygon Boolean operations](polybool2d.md) — done (Wave 3).
- [Non-rigid ICP registration](nonrigid.md) — done (Wave 3).
- [2D convex hull](hull2d.md) — done (Wave 3).
- [2D convex Minkowski sum](minkowski2d.md) — done (Wave 3).
- [Polygon triangulation with holes](polytri_holes.md) — done (Wave 3).
- [Polyline simplification](simplify.md) — done (Wave 3).
- [MLS surface projection](mls.md) — done (Wave 3).
- [SDF ray casting](sdf_raycast.md) — done (Wave 3).
- [Dijkstra geodesic distances](geodesic_dijkstra.md) — done (Wave 3).
- [Angle-based flattening](abf.md) — done (Wave 3).
- [Mean-value cage deformation](meancoords.md) — done (Wave 3).
- [Bilateral mesh smoothing](mesh_bilateral.md) — done (Wave 3).
- [Catmull-Rom spline resampling](spline.md) — done (Wave 3).
- [Surface mesh edge-flip quality improvement](quality_flip.md) — done (Wave 3).
- [Isotropic edge-split refinement](remesh_split.md) — done (Wave 3).
- [Isotropic edge-collapse coarsening](remesh_collapse.md) — done (Wave 3).
- [2D Delaunay refinement](refine2d.md) — done (Wave 3).
- [Surface Voronoi segmentation](surface_voronoi.md) — done (Wave 3).
- [Point-to-plane ICP registration](icp_plane.md) — done (Wave 3).
- [Mesh contour lines](mesh_contour.md) — done (Wave 3).
- [Dihedral-angle region growing](region_growing.md) — done (Wave 3).
- [Point-cloud local curvature](point_curvature.md) — done (Wave 3).
- [Heat kernel signature](hks.md) — done (Wave 3).
- [2D convex polygon decomposition](convex_decompose2d.md) — done (Wave 3).
- [Ritter bounding sphere](bounding_sphere.md) — done (Wave 3).
- [Minimum enclosing sphere](miniball.md) — done (Wave 3).
- [Mesh mass properties](mass_properties.md) — done (Wave 3).
- [Largest inscribed circle](inscribed_circle.md) — done (Wave 3).
- [Minimum-area bounding rectangle](min_area_rect.md) — done (Wave 3).
- [Polyline-mesh intersection](polyline_intersection.md) — done (Wave 3).
- [Mesh ambient occlusion](ambient_occlusion.md) — done (Wave 3).
- [Spectral mesh segmentation](spectral_segment.md) — done (Wave 3).
- [Mesh-mesh minimum distance](mesh_distance.md) — done (Wave 3).
- [Reflection symmetry detection](symmetry.md) — done (Wave 3).
- [Point-cloud k-means](kmeans.md) — done (Wave 3).
- [DBSCAN clustering](dbscan.md) — done (Wave 3).
- [Arc-length resampling](resample.md) — done (Wave 3).
- [2D polygon moments](polygon_moments.md) — done (Wave 3).
- [Convex polygon clipping](polygon_clip.md) — done (Wave 3).
- [Exact 2D minimum enclosing circle](min_circle2d.md) — done (Wave 3).
- [Geodesic shortest-path polyline](geodesic_path.md) — done (Wave 3).
- [Geodesic ball extraction](geodesic_ball.md) — done (Wave 3).
- [Rigid point-set fitting](rigid_fit.md) — done (Wave 3).
- [2D least-squares circle fitting](circle_fit.md) — done (Wave 3).
- [3D least-squares sphere fitting](sphere_fit.md) — done (Wave 3).
- [3D least-squares plane fitting](plane_fit.md) — done (Wave 3).
- [Chaikin polyline subdivision](chaikin.md) — done (Wave 3).
- [2D ellipse fitting](ellipse_fit.md) — done (Wave 3).
- [Least-squares cylinder fitting](cylinder_fit.md) — done (Wave 3).
- [3D segment-segment distance](segment_distance.md) — done (Wave 3).
- [Ray/plane intersection](ray_plane.md) — done (Wave 3).
- [Triangle-triangle intersection](tri_tri.md) — done (Wave 3).
- [Polygon Boolean with holes](polybool_holes.md) — done (Wave 3).
- [Geodesic diameter](geodesic_diameter.md) — done (Wave 3).
- [Tetrahedral cage deformation](cage3d.md) — done (Wave 3).
- [3D volume Poisson-disk sampling](poisson3d.md) — done (Wave 3).
- [2D point-in-polygon](point_in_polygon.md) — done (Wave 3).
- [Cotangent Laplacian](laplacian.md) — done (Wave 3).
- [Quaternion operations](quaternion.md) — done (Wave 3).
- [Mesh plane clipping](mesh_clip.md) — done (Wave 3).
- [2D segment intersection](segment2d.md) — done (Wave 3).
- [Minimum spanning tree](mst.md) — done (Wave 3).
- [Principal component analysis](pca.md) — done (Wave 3).
- [Sphere–sphere intersection](sphere_intersection.md) — done (Wave 3).
- [Discrete Fréchet distance](frechet.md) — done (Wave 3).
- [2D polygon Poisson-disk sampling](poisson2d.md) — done (Wave 3).
- [Iterative Laplacian curve smoothing](curve_smooth.md) — done (Wave 3).
- [2D harmonic cage coordinates](cage2d_harmonic.md) — done (Wave 3).
- [2D polygon rasterization](raster2d.md) — done (Wave 3).
- [Visibility-graph shortest path](visibility.md) — done (Wave 3).
- [Broad-phase AABB overlap pairs](mesh_overlap.md) — done (Wave 3).
- [Solid angle and winding number](solid_angle.md) — done (Wave 3).
- [Visvalingam-Whyatt simplification](vw.md) — done (Wave 3).
- [Bezier evaluation and subdivision](bezier.md) — done (Wave 3).
- [Ray–sphere intersection](ray_sphere.md) — done (Wave 3).
- [Plane–plane intersection](plane_plane.md) — done (Wave 3).
- [Sphere–AABB overlap](sphere_aabb.md) — done (Wave 3).
- [2D circle–circle intersection](circle_circle.md) — done (Wave 3).
- [Segment–triangle intersection](segment_triangle.md) — done (Wave 3).
- [Ray–cylinder intersection](ray_cylinder.md) — done (Wave 3).
- [Capsule collision](capsule.md) — done (Wave 3).
- [Sphere–triangle overlap](sphere_triangle.md) — done (Wave 3).
- [Triangle–AABB overlap](tri_aabb.md) — done (Wave 3).
- [Convex polygon diameter](polygon_diameter.md) — done (Wave 3).
- [Convex polygon minimum width](polygon_width.md) — done (Wave 3).
- [AABB–AABB distance](aabb_distance.md) — done (Wave 3).
- [Ray–OBB intersection](ray_obb.md) — done (Wave 3).
- [Triangle–triangle distance](tri_tri_distance.md) — done (Wave 3).
- [Plane–sphere intersection](plane_sphere.md) — done (Wave 3).
- [2D circle–line intersection](circle_line.md) — done (Wave 3).
- [Ray–tetrahedron intersection](ray_tetra.md) — done (Wave 3).
- [Closest point on a polyline](polyline_closest.md) — done (Wave 3).
- [Capsule–AABB overlap](capsule_aabb.md) — done (Wave 3).
- [2D polygon–polygon distance](polygon_distance2d.md) — done (Wave 3).
- [Triangle–tetrahedron intersection](tri_tetra.md) — done (Wave 3).
- [Sphere–capsule overlap](sphere_capsule.md) — done (Wave 3).
- [Tetrahedron–tetrahedron intersection](tetra_tetra.md) — done (Wave 3).
- [Sphere–OBB overlap](sphere_obb.md) — done (Wave 3).
- [2D circle–AABB overlap](circle_aabb.md) — done (Wave 3).
- [Capsule–OBB overlap](capsule_obb.md) — done (Wave 3).
- [Ray–cone intersection](ray_cone.md) — done (Wave 3).
- [Triangle–OBB overlap](tri_obb.md) — done (Wave 3).
- [AABB–OBB overlap](aabb_obb.md) — done (Wave 3).
- [OBB–OBB overlap](obb_obb.md) — done (Wave 3).
- [Ellipse point queries](ellipse.md) — done (Wave 3).
- [Capsule–triangle overlap](capsule_triangle.md) — done (Wave 3).
- [Closest point on an ellipse](ellipse_closest.md) — done (Wave 3).
- [Ellipse area and perimeter](ellipse_measures.md) — done (Wave 3).
- [2D circle–triangle intersection](circle_triangle.md) — done (Wave 3).
- [2D line–line intersection](line2d.md) — done (Wave 3).
- [Polyline arc length](arc_length.md) — done (Wave 3).
- [Circle intersection and union area](circle_area.md) — done (Wave 3).
- [Convex polygon intersection area](polygon_intersection_area.md) — done (Wave 3).
- [Circle–ellipse overlap](circle_ellipse.md) — done (Wave 3).
- [MMP continuous geodesic](geodesic_mmp.md) — done (Wave 3).
- [Anisotropic metric tensors](metric.md) — done (Wave 3).
- Future: anisotropic metric-tensor remeshing (restricted Voronoi diagrams).
