add_executable(test_linalg test_linalg.cpp)
target_link_libraries(test_linalg PRIVATE geomlib)
add_test(NAME test_linalg COMMAND test_linalg)

add_executable(test_linalg_extended test_linalg_extended.cpp)
target_link_libraries(test_linalg_extended PRIVATE geomlib)
add_test(NAME test_linalg_extended COMMAND test_linalg_extended)

add_executable(test_solver test_solver.cpp)
target_link_libraries(test_solver PRIVATE geomlib)
add_test(NAME test_solver COMMAND test_solver)

add_executable(test_mesh test_mesh.cpp)
target_link_libraries(test_mesh PRIVATE geomlib)
add_test(NAME test_mesh COMMAND test_mesh)

add_executable(test_mesh_extended test_mesh_extended.cpp)
target_link_libraries(test_mesh_extended PRIVATE geomlib)
add_test(NAME test_mesh_extended COMMAND test_mesh_extended)

add_executable(test_random test_random.cpp)
target_link_libraries(test_random PRIVATE geomlib)
add_test(NAME test_random COMMAND test_random)

add_executable(test_montecarlo test_montecarlo.cpp)
target_link_libraries(test_montecarlo PRIVATE geomlib)
add_test(NAME test_montecarlo COMMAND test_montecarlo)

add_executable(test_bvh test_bvh.cpp)
target_link_libraries(test_bvh PRIVATE geomlib)
add_test(NAME test_bvh COMMAND test_bvh)

add_executable(test_bvh_extended test_bvh_extended.cpp)
target_link_libraries(test_bvh_extended PRIVATE geomlib)
add_test(NAME test_bvh_extended COMMAND test_bvh_extended)

add_executable(test_integration test_integration.cpp)
target_link_libraries(test_integration PRIVATE geomlib)
add_test(NAME test_integration COMMAND test_integration)

add_executable(test_collision test_collision.cpp)
target_link_libraries(test_collision PRIVATE geomlib)
add_test(NAME test_collision COMMAND test_collision)

add_executable(test_advanced_features test_advanced_features.cpp)
target_link_libraries(test_advanced_features PRIVATE geomlib)
add_test(NAME test_advanced_features COMMAND test_advanced_features)

add_executable(test_optimize test_optimize.cpp)
target_link_libraries(test_optimize PRIVATE geomlib)
add_test(NAME test_optimize COMMAND test_optimize)

add_executable(test_cfd test_cfd.cpp)
target_link_libraries(test_cfd PRIVATE geomlib)
add_test(NAME test_cfd COMMAND test_cfd)

add_executable(test_sdf test_sdf.cpp)
target_link_libraries(test_sdf PRIVATE geomlib)
add_test(NAME test_sdf COMMAND test_sdf)

add_executable(test_smoothing test_smoothing.cpp)
target_link_libraries(test_smoothing PRIVATE geomlib)
add_test(NAME test_smoothing COMMAND test_smoothing)

add_executable(test_subdivision test_subdivision.cpp)
target_link_libraries(test_subdivision PRIVATE geomlib)
add_test(NAME test_subdivision COMMAND test_subdivision)

add_executable(test_geodesic test_geodesic.cpp)
target_link_libraries(test_geodesic PRIVATE geomlib)
add_test(NAME test_geodesic COMMAND test_geodesic)

add_executable(test_reconstruction test_reconstruction.cpp)
target_link_libraries(test_reconstruction PRIVATE geomlib)
add_test(NAME test_reconstruction COMMAND test_reconstruction)

add_executable(test_decimation test_decimation.cpp)
target_link_libraries(test_decimation PRIVATE geomlib)
add_test(NAME test_decimation COMMAND test_decimation)

add_executable(test_curvature test_curvature.cpp)
target_link_libraries(test_curvature PRIVATE geomlib)
add_test(NAME test_curvature COMMAND test_curvature)

add_executable(test_atlas test_atlas.cpp)
target_link_libraries(test_atlas PRIVATE geomlib)
add_test(NAME test_atlas COMMAND test_atlas)

add_executable(test_overdraw test_overdraw.cpp)
target_link_libraries(test_overdraw PRIVATE geomlib)
add_test(NAME test_overdraw COMMAND test_overdraw)

add_executable(test_predicates test_predicates.cpp)
target_link_libraries(test_predicates PRIVATE geomlib)
add_test(NAME test_predicates COMMAND test_predicates)

add_executable(test_sparse test_sparse.cpp)
target_link_libraries(test_sparse PRIVATE geomlib)
add_test(NAME test_sparse COMMAND test_sparse)

add_executable(test_io test_io.cpp)
target_link_libraries(test_io PRIVATE geomlib)
add_test(NAME test_io COMMAND test_io)

add_executable(test_distance test_distance.cpp)
target_link_libraries(test_distance PRIVATE geomlib)
add_test(NAME test_distance COMMAND test_distance)

add_executable(test_parameterization test_parameterization.cpp)
target_link_libraries(test_parameterization PRIVATE geomlib)
add_test(NAME test_parameterization COMMAND test_parameterization)

add_executable(test_repair test_repair.cpp)
target_link_libraries(test_repair PRIVATE geomlib)
add_test(NAME test_repair COMMAND test_repair)

add_executable(test_boolean test_boolean.cpp)
target_link_libraries(test_boolean PRIVATE geomlib)
add_test(NAME test_boolean COMMAND test_boolean)

add_executable(test_registration test_registration.cpp)
target_link_libraries(test_registration PRIVATE geomlib)
add_test(NAME test_registration COMMAND test_registration)

add_executable(test_bluenoise test_bluenoise.cpp)
target_link_libraries(test_bluenoise PRIVATE geomlib)
add_test(NAME test_bluenoise COMMAND test_bluenoise)

add_executable(test_poisson test_poisson.cpp)
target_link_libraries(test_poisson PRIVATE geomlib)
add_test(NAME test_poisson COMMAND test_poisson)

add_executable(test_acvd test_acvd.cpp)
target_link_libraries(test_acvd PRIVATE geomlib)
add_test(NAME test_acvd COMMAND test_acvd)

add_executable(test_deform test_deform.cpp)
target_link_libraries(test_deform PRIVATE geomlib)
add_test(NAME test_deform COMMAND test_deform)

add_executable(test_harmonics test_harmonics.cpp)
target_link_libraries(test_harmonics PRIVATE geomlib)
add_test(NAME test_harmonics COMMAND test_harmonics)

add_executable(test_shape_detection test_shape_detection.cpp)
target_link_libraries(test_shape_detection PRIVATE geomlib)
add_test(NAME test_shape_detection COMMAND test_shape_detection)

add_executable(test_quadremesh test_quadremesh.cpp)
target_link_libraries(test_quadremesh PRIVATE geomlib)
add_test(NAME test_quadremesh COMMAND test_quadremesh)

add_executable(test_anisotropic test_anisotropic.cpp)
target_link_libraries(test_anisotropic PRIVATE geomlib)
add_test(NAME test_anisotropic COMMAND test_anisotropic)

add_executable(test_marching_tets test_marching_tets.cpp)
target_link_libraries(test_marching_tets PRIVATE geomlib)
add_test(NAME test_marching_tets COMMAND test_marching_tets)

add_executable(test_dual_contouring test_dual_contouring.cpp)
target_link_libraries(test_dual_contouring PRIVATE geomlib)
add_test(NAME test_dual_contouring COMMAND test_dual_contouring)

add_executable(test_saliency test_saliency.cpp)
target_link_libraries(test_saliency PRIVATE geomlib)
add_test(NAME test_saliency COMMAND test_saliency)

add_executable(test_skeleton test_skeleton.cpp)
target_link_libraries(test_skeleton PRIVATE geomlib)
add_test(NAME test_skeleton COMMAND test_skeleton)

add_executable(test_straight_skeleton test_straight_skeleton.cpp)
target_link_libraries(test_straight_skeleton PRIVATE geomlib)
add_test(NAME test_straight_skeleton COMMAND test_straight_skeleton)

add_executable(test_polygon_offset test_polygon_offset.cpp)
target_link_libraries(test_polygon_offset PRIVATE geomlib)
add_test(NAME test_polygon_offset COMMAND test_polygon_offset)

add_executable(test_vsa test_vsa.cpp)
target_link_libraries(test_vsa PRIVATE geomlib)
add_test(NAME test_vsa COMMAND test_vsa)

add_executable(test_pointset test_pointset.cpp)
target_link_libraries(test_pointset PRIVATE geomlib)
add_test(NAME test_pointset COMMAND test_pointset)

add_executable(test_alpha_wrap test_alpha_wrap.cpp)
target_link_libraries(test_alpha_wrap PRIVATE geomlib)
add_test(NAME test_alpha_wrap COMMAND test_alpha_wrap)

add_executable(test_slicer test_slicer.cpp)
target_link_libraries(test_slicer PRIVATE geomlib)
add_test(NAME test_slicer COMMAND test_slicer)

add_executable(test_convex_decomposition test_convex_decomposition.cpp)
target_link_libraries(test_convex_decomposition PRIVATE geomlib)
add_test(NAME test_convex_decomposition COMMAND test_convex_decomposition)

add_executable(test_fairing test_fairing.cpp)
target_link_libraries(test_fairing PRIVATE geomlib)
add_test(NAME test_fairing COMMAND test_fairing)

add_executable(test_fast_winding test_fast_winding.cpp)
target_link_libraries(test_fast_winding PRIVATE geomlib)
add_test(NAME test_fast_winding COMMAND test_fast_winding)

add_executable(test_partition test_partition.cpp)
target_link_libraries(test_partition PRIVATE geomlib)
add_test(NAME test_partition COMMAND test_partition)

add_executable(test_clustering test_clustering.cpp)
target_link_libraries(test_clustering PRIVATE geomlib)
add_test(NAME test_clustering COMMAND test_clustering)

add_executable(test_marching_squares test_marching_squares.cpp)
target_link_libraries(test_marching_squares PRIVATE geomlib)
add_test(NAME test_marching_squares COMMAND test_marching_squares)

add_executable(test_delaunay2d test_delaunay2d.cpp)
target_link_libraries(test_delaunay2d PRIVATE geomlib)
add_test(NAME test_delaunay2d COMMAND test_delaunay2d)

add_executable(test_geodesic_fmm test_geodesic_fmm.cpp)
target_link_libraries(test_geodesic_fmm PRIVATE geomlib)
add_test(NAME test_geodesic_fmm COMMAND test_geodesic_fmm)

add_executable(test_geodesic_exact test_geodesic_exact.cpp)
target_link_libraries(test_geodesic_exact PRIVATE geomlib)
add_test(NAME test_geodesic_exact COMMAND test_geodesic_exact)

add_executable(test_skinning test_skinning.cpp)
target_link_libraries(test_skinning PRIVATE geomlib)
add_test(NAME test_skinning COMMAND test_skinning)

add_executable(test_dualquat_skinning test_dualquat_skinning.cpp)
target_link_libraries(test_dualquat_skinning PRIVATE geomlib)
add_test(NAME test_dualquat_skinning COMMAND test_dualquat_skinning)

add_executable(test_cross_field test_cross_field.cpp)
target_link_libraries(test_cross_field PRIVATE geomlib)
add_test(NAME test_cross_field COMMAND test_cross_field)

add_executable(test_shape_diameter test_shape_diameter.cpp)
target_link_libraries(test_shape_diameter PRIVATE geomlib)
add_test(NAME test_shape_diameter COMMAND test_shape_diameter)

add_executable(test_kdtree test_kdtree.cpp)
target_link_libraries(test_kdtree PRIVATE geomlib)
add_test(NAME test_kdtree COMMAND test_kdtree)

add_executable(test_bvh_refit test_bvh_refit.cpp)
target_link_libraries(test_bvh_refit PRIVATE geomlib)
add_test(NAME test_bvh_refit COMMAND test_bvh_refit)

add_executable(test_octree test_octree.cpp)
target_link_libraries(test_octree PRIVATE geomlib)
add_test(NAME test_octree COMMAND test_octree)

add_executable(test_sparse_direct test_sparse_direct.cpp)
target_link_libraries(test_sparse_direct PRIVATE geomlib)
add_test(NAME test_sparse_direct COMMAND test_sparse_direct)

add_executable(test_hexdominant test_hexdominant.cpp)
target_link_libraries(test_hexdominant PRIVATE geomlib)
add_test(NAME test_hexdominant COMMAND test_hexdominant)

add_executable(test_rbf_reconstruct test_rbf_reconstruct.cpp)
target_link_libraries(test_rbf_reconstruct PRIVATE geomlib)
add_test(NAME test_rbf_reconstruct COMMAND test_rbf_reconstruct)

add_executable(test_fluid test_fluid.cpp)
target_link_libraries(test_fluid PRIVATE geomlib)
add_test(NAME test_fluid COMMAND test_fluid)

add_executable(test_voronoi2d test_voronoi2d.cpp)
target_link_libraries(test_voronoi2d PRIVATE geomlib)
add_test(NAME test_voronoi2d COMMAND test_voronoi2d)

add_executable(test_voxelize test_voxelize.cpp)
target_link_libraries(test_voxelize PRIVATE geomlib)
add_test(NAME test_voxelize COMMAND test_voxelize)

add_executable(test_spatial_sort test_spatial_sort.cpp)
target_link_libraries(test_spatial_sort PRIVATE geomlib)
add_test(NAME test_spatial_sort COMMAND test_spatial_sort)

add_executable(test_intersection test_intersection.cpp)
target_link_libraries(test_intersection PRIVATE geomlib)
add_test(NAME test_intersection COMMAND test_intersection)

add_executable(test_cdt2d test_cdt2d.cpp)
target_link_libraries(test_cdt2d PRIVATE geomlib)
add_test(NAME test_cdt2d COMMAND test_cdt2d)

add_executable(test_offset test_offset.cpp)
target_link_libraries(test_offset PRIVATE geomlib)
add_test(NAME test_offset COMMAND test_offset)

add_executable(test_fps test_fps.cpp)
target_link_libraries(test_fps PRIVATE geomlib)
add_test(NAME test_fps COMMAND test_fps)

add_executable(test_medial2d test_medial2d.cpp)
target_link_libraries(test_medial2d PRIVATE geomlib)
add_test(NAME test_medial2d COMMAND test_medial2d)

add_executable(test_polybool2d test_polybool2d.cpp)
target_link_libraries(test_polybool2d PRIVATE geomlib)
add_test(NAME test_polybool2d COMMAND test_polybool2d)

add_executable(test_nonrigid test_nonrigid.cpp)
target_link_libraries(test_nonrigid PRIVATE geomlib)
add_test(NAME test_nonrigid COMMAND test_nonrigid)

add_executable(test_hull2d test_hull2d.cpp)
target_link_libraries(test_hull2d PRIVATE geomlib)
add_test(NAME test_hull2d COMMAND test_hull2d)

add_executable(test_minkowski2d test_minkowski2d.cpp)
target_link_libraries(test_minkowski2d PRIVATE geomlib)
add_test(NAME test_minkowski2d COMMAND test_minkowski2d)

add_executable(test_polytri_holes test_polytri_holes.cpp)
target_link_libraries(test_polytri_holes PRIVATE geomlib)
add_test(NAME test_polytri_holes COMMAND test_polytri_holes)

add_executable(test_simplify test_simplify.cpp)
target_link_libraries(test_simplify PRIVATE geomlib)
add_test(NAME test_simplify COMMAND test_simplify)

add_executable(test_mls test_mls.cpp)
target_link_libraries(test_mls PRIVATE geomlib)
add_test(NAME test_mls COMMAND test_mls)

add_executable(test_sdf_raycast test_sdf_raycast.cpp)
target_link_libraries(test_sdf_raycast PRIVATE geomlib)
add_test(NAME test_sdf_raycast COMMAND test_sdf_raycast)

add_executable(test_geodesic_dijkstra test_geodesic_dijkstra.cpp)
target_link_libraries(test_geodesic_dijkstra PRIVATE geomlib)
add_test(NAME test_geodesic_dijkstra COMMAND test_geodesic_dijkstra)

add_executable(test_abf test_abf.cpp)
target_link_libraries(test_abf PRIVATE geomlib)
add_test(NAME test_abf COMMAND test_abf)

add_executable(test_meancoords test_meancoords.cpp)
target_link_libraries(test_meancoords PRIVATE geomlib)
add_test(NAME test_meancoords COMMAND test_meancoords)

add_executable(test_mesh_bilateral test_mesh_bilateral.cpp)
target_link_libraries(test_mesh_bilateral PRIVATE geomlib)
add_test(NAME test_mesh_bilateral COMMAND test_mesh_bilateral)

add_executable(test_spline test_spline.cpp)
target_link_libraries(test_spline PRIVATE geomlib)
add_test(NAME test_spline COMMAND test_spline)

add_executable(test_quality_flip test_quality_flip.cpp)
target_link_libraries(test_quality_flip PRIVATE geomlib)
add_test(NAME test_quality_flip COMMAND test_quality_flip)

add_executable(test_remesh_split test_remesh_split.cpp)
target_link_libraries(test_remesh_split PRIVATE geomlib)
add_test(NAME test_remesh_split COMMAND test_remesh_split)

add_executable(test_remesh_collapse test_remesh_collapse.cpp)
target_link_libraries(test_remesh_collapse PRIVATE geomlib)
add_test(NAME test_remesh_collapse COMMAND test_remesh_collapse)

add_executable(test_refine2d test_refine2d.cpp)
target_link_libraries(test_refine2d PRIVATE geomlib)
add_test(NAME test_refine2d COMMAND test_refine2d)

add_executable(test_surface_voronoi test_surface_voronoi.cpp)
target_link_libraries(test_surface_voronoi PRIVATE geomlib)
add_test(NAME test_surface_voronoi COMMAND test_surface_voronoi)

add_executable(test_icp_plane test_icp_plane.cpp)
target_link_libraries(test_icp_plane PRIVATE geomlib)
add_test(NAME test_icp_plane COMMAND test_icp_plane)

add_executable(test_mesh_contour test_mesh_contour.cpp)
target_link_libraries(test_mesh_contour PRIVATE geomlib)
add_test(NAME test_mesh_contour COMMAND test_mesh_contour)

add_executable(test_region_growing test_region_growing.cpp)
target_link_libraries(test_region_growing PRIVATE geomlib)
add_test(NAME test_region_growing COMMAND test_region_growing)

add_executable(test_point_curvature test_point_curvature.cpp)
target_link_libraries(test_point_curvature PRIVATE geomlib)
add_test(NAME test_point_curvature COMMAND test_point_curvature)

add_executable(test_hks test_hks.cpp)
target_link_libraries(test_hks PRIVATE geomlib)
add_test(NAME test_hks COMMAND test_hks)

add_executable(test_convex_decompose2d test_convex_decompose2d.cpp)
target_link_libraries(test_convex_decompose2d PRIVATE geomlib)
add_test(NAME test_convex_decompose2d COMMAND test_convex_decompose2d)

add_executable(test_bounding_sphere test_bounding_sphere.cpp)
target_link_libraries(test_bounding_sphere PRIVATE geomlib)
add_test(NAME test_bounding_sphere COMMAND test_bounding_sphere)

add_executable(test_miniball test_miniball.cpp)
target_link_libraries(test_miniball PRIVATE geomlib)
add_test(NAME test_miniball COMMAND test_miniball)

add_executable(test_mass_properties test_mass_properties.cpp)
target_link_libraries(test_mass_properties PRIVATE geomlib)
add_test(NAME test_mass_properties COMMAND test_mass_properties)

add_executable(test_inscribed_circle test_inscribed_circle.cpp)
target_link_libraries(test_inscribed_circle PRIVATE geomlib)
add_test(NAME test_inscribed_circle COMMAND test_inscribed_circle)

add_executable(test_min_area_rect test_min_area_rect.cpp)
target_link_libraries(test_min_area_rect PRIVATE geomlib)
add_test(NAME test_min_area_rect COMMAND test_min_area_rect)

add_executable(test_polyline_intersection test_polyline_intersection.cpp)
target_link_libraries(test_polyline_intersection PRIVATE geomlib)
add_test(NAME test_polyline_intersection COMMAND test_polyline_intersection)

add_executable(test_ambient_occlusion test_ambient_occlusion.cpp)
target_link_libraries(test_ambient_occlusion PRIVATE geomlib)
add_test(NAME test_ambient_occlusion COMMAND test_ambient_occlusion)

add_executable(test_spectral_segment test_spectral_segment.cpp)
target_link_libraries(test_spectral_segment PRIVATE geomlib)
add_test(NAME test_spectral_segment COMMAND test_spectral_segment)

add_executable(test_mesh_distance test_mesh_distance.cpp)
target_link_libraries(test_mesh_distance PRIVATE geomlib)
add_test(NAME test_mesh_distance COMMAND test_mesh_distance)

add_executable(test_symmetry test_symmetry.cpp)
target_link_libraries(test_symmetry PRIVATE geomlib)
add_test(NAME test_symmetry COMMAND test_symmetry)

add_executable(test_kmeans test_kmeans.cpp)
target_link_libraries(test_kmeans PRIVATE geomlib)
add_test(NAME test_kmeans COMMAND test_kmeans)

add_executable(test_dbscan test_dbscan.cpp)
target_link_libraries(test_dbscan PRIVATE geomlib)
add_test(NAME test_dbscan COMMAND test_dbscan)

add_executable(test_resample test_resample.cpp)
target_link_libraries(test_resample PRIVATE geomlib)
add_test(NAME test_resample COMMAND test_resample)

add_executable(test_polygon_moments test_polygon_moments.cpp)
target_link_libraries(test_polygon_moments PRIVATE geomlib)
add_test(NAME test_polygon_moments COMMAND test_polygon_moments)

add_executable(test_polygon_clip test_polygon_clip.cpp)
target_link_libraries(test_polygon_clip PRIVATE geomlib)
add_test(NAME test_polygon_clip COMMAND test_polygon_clip)

add_executable(test_min_circle2d test_min_circle2d.cpp)
target_link_libraries(test_min_circle2d PRIVATE geomlib)
add_test(NAME test_min_circle2d COMMAND test_min_circle2d)

add_executable(test_geodesic_path test_geodesic_path.cpp)
target_link_libraries(test_geodesic_path PRIVATE geomlib)
add_test(NAME test_geodesic_path COMMAND test_geodesic_path)

add_executable(test_geodesic_ball test_geodesic_ball.cpp)
target_link_libraries(test_geodesic_ball PRIVATE geomlib)
add_test(NAME test_geodesic_ball COMMAND test_geodesic_ball)

add_executable(test_rigid_fit test_rigid_fit.cpp)
target_link_libraries(test_rigid_fit PRIVATE geomlib)
add_test(NAME test_rigid_fit COMMAND test_rigid_fit)

add_executable(test_circle_fit test_circle_fit.cpp)
target_link_libraries(test_circle_fit PRIVATE geomlib)
add_test(NAME test_circle_fit COMMAND test_circle_fit)

add_executable(test_sphere_fit test_sphere_fit.cpp)
target_link_libraries(test_sphere_fit PRIVATE geomlib)
add_test(NAME test_sphere_fit COMMAND test_sphere_fit)

add_executable(test_plane_fit test_plane_fit.cpp)
target_link_libraries(test_plane_fit PRIVATE geomlib)
add_test(NAME test_plane_fit COMMAND test_plane_fit)

add_executable(test_chaikin test_chaikin.cpp)
target_link_libraries(test_chaikin PRIVATE geomlib)
add_test(NAME test_chaikin COMMAND test_chaikin)

add_executable(test_ellipse_fit test_ellipse_fit.cpp)
target_link_libraries(test_ellipse_fit PRIVATE geomlib)
add_test(NAME test_ellipse_fit COMMAND test_ellipse_fit)

add_executable(test_cylinder_fit test_cylinder_fit.cpp)
target_link_libraries(test_cylinder_fit PRIVATE geomlib)
add_test(NAME test_cylinder_fit COMMAND test_cylinder_fit)

add_executable(test_segment_distance test_segment_distance.cpp)
target_link_libraries(test_segment_distance PRIVATE geomlib)
add_test(NAME test_segment_distance COMMAND test_segment_distance)

add_executable(test_ray_plane test_ray_plane.cpp)
target_link_libraries(test_ray_plane PRIVATE geomlib)
add_test(NAME test_ray_plane COMMAND test_ray_plane)

add_executable(test_tri_tri test_tri_tri.cpp)
target_link_libraries(test_tri_tri PRIVATE geomlib)
add_test(NAME test_tri_tri COMMAND test_tri_tri)

add_executable(test_polybool_holes test_polybool_holes.cpp)
target_link_libraries(test_polybool_holes PRIVATE geomlib)
add_test(NAME test_polybool_holes COMMAND test_polybool_holes)

add_executable(test_geodesic_diameter test_geodesic_diameter.cpp)
target_link_libraries(test_geodesic_diameter PRIVATE geomlib)
add_test(NAME test_geodesic_diameter COMMAND test_geodesic_diameter)

add_executable(test_cage3d test_cage3d.cpp)
target_link_libraries(test_cage3d PRIVATE geomlib)
add_test(NAME test_cage3d COMMAND test_cage3d)

add_executable(test_poisson3d test_poisson3d.cpp)
target_link_libraries(test_poisson3d PRIVATE geomlib)
add_test(NAME test_poisson3d COMMAND test_poisson3d)

add_executable(test_point_in_polygon test_point_in_polygon.cpp)
target_link_libraries(test_point_in_polygon PRIVATE geomlib)
add_test(NAME test_point_in_polygon COMMAND test_point_in_polygon)

add_executable(test_laplacian test_laplacian.cpp)
target_link_libraries(test_laplacian PRIVATE geomlib)
add_test(NAME test_laplacian COMMAND test_laplacian)

add_executable(test_quaternion test_quaternion.cpp)
target_link_libraries(test_quaternion PRIVATE geomlib)
add_test(NAME test_quaternion COMMAND test_quaternion)

add_executable(test_mesh_clip test_mesh_clip.cpp)
target_link_libraries(test_mesh_clip PRIVATE geomlib)
add_test(NAME test_mesh_clip COMMAND test_mesh_clip)

add_executable(test_segment2d test_segment2d.cpp)
target_link_libraries(test_segment2d PRIVATE geomlib)
add_test(NAME test_segment2d COMMAND test_segment2d)

add_executable(test_mst test_mst.cpp)
target_link_libraries(test_mst PRIVATE geomlib)
add_test(NAME test_mst COMMAND test_mst)

add_executable(test_pca test_pca.cpp)
target_link_libraries(test_pca PRIVATE geomlib)
add_test(NAME test_pca COMMAND test_pca)

add_executable(test_sphere_intersection test_sphere_intersection.cpp)
target_link_libraries(test_sphere_intersection PRIVATE geomlib)
add_test(NAME test_sphere_intersection COMMAND test_sphere_intersection)

add_executable(test_frechet test_frechet.cpp)
target_link_libraries(test_frechet PRIVATE geomlib)
add_test(NAME test_frechet COMMAND test_frechet)

add_executable(test_poisson2d test_poisson2d.cpp)
target_link_libraries(test_poisson2d PRIVATE geomlib)
add_test(NAME test_poisson2d COMMAND test_poisson2d)

add_executable(test_curve_smooth test_curve_smooth.cpp)
target_link_libraries(test_curve_smooth PRIVATE geomlib)
add_test(NAME test_curve_smooth COMMAND test_curve_smooth)

add_executable(test_cage2d_harmonic test_cage2d_harmonic.cpp)
target_link_libraries(test_cage2d_harmonic PRIVATE geomlib)
add_test(NAME test_cage2d_harmonic COMMAND test_cage2d_harmonic)

add_executable(test_raster2d test_raster2d.cpp)
target_link_libraries(test_raster2d PRIVATE geomlib)
add_test(NAME test_raster2d COMMAND test_raster2d)

add_executable(test_visibility test_visibility.cpp)
target_link_libraries(test_visibility PRIVATE geomlib)
add_test(NAME test_visibility COMMAND test_visibility)

add_executable(test_mesh_overlap test_mesh_overlap.cpp)
target_link_libraries(test_mesh_overlap PRIVATE geomlib)
add_test(NAME test_mesh_overlap COMMAND test_mesh_overlap)

add_executable(test_solid_angle test_solid_angle.cpp)
target_link_libraries(test_solid_angle PRIVATE geomlib)
add_test(NAME test_solid_angle COMMAND test_solid_angle)

add_executable(test_vw test_vw.cpp)
target_link_libraries(test_vw PRIVATE geomlib)
add_test(NAME test_vw COMMAND test_vw)

add_executable(test_bezier test_bezier.cpp)
target_link_libraries(test_bezier PRIVATE geomlib)
add_test(NAME test_bezier COMMAND test_bezier)

add_executable(test_ray_sphere test_ray_sphere.cpp)
target_link_libraries(test_ray_sphere PRIVATE geomlib)
add_test(NAME test_ray_sphere COMMAND test_ray_sphere)

add_executable(test_plane_plane test_plane_plane.cpp)
target_link_libraries(test_plane_plane PRIVATE geomlib)
add_test(NAME test_plane_plane COMMAND test_plane_plane)

add_executable(test_sphere_aabb test_sphere_aabb.cpp)
target_link_libraries(test_sphere_aabb PRIVATE geomlib)
add_test(NAME test_sphere_aabb COMMAND test_sphere_aabb)

add_executable(test_circle_circle test_circle_circle.cpp)
target_link_libraries(test_circle_circle PRIVATE geomlib)
add_test(NAME test_circle_circle COMMAND test_circle_circle)

add_executable(test_segment_triangle test_segment_triangle.cpp)
target_link_libraries(test_segment_triangle PRIVATE geomlib)
add_test(NAME test_segment_triangle COMMAND test_segment_triangle)

add_executable(test_ray_cylinder test_ray_cylinder.cpp)
target_link_libraries(test_ray_cylinder PRIVATE geomlib)
add_test(NAME test_ray_cylinder COMMAND test_ray_cylinder)

add_executable(test_capsule test_capsule.cpp)
target_link_libraries(test_capsule PRIVATE geomlib)
add_test(NAME test_capsule COMMAND test_capsule)

add_executable(test_sphere_triangle test_sphere_triangle.cpp)
target_link_libraries(test_sphere_triangle PRIVATE geomlib)
add_test(NAME test_sphere_triangle COMMAND test_sphere_triangle)

add_executable(test_tri_aabb test_tri_aabb.cpp)
target_link_libraries(test_tri_aabb PRIVATE geomlib)
add_test(NAME test_tri_aabb COMMAND test_tri_aabb)

add_executable(test_polygon_diameter test_polygon_diameter.cpp)
target_link_libraries(test_polygon_diameter PRIVATE geomlib)
add_test(NAME test_polygon_diameter COMMAND test_polygon_diameter)

add_executable(test_polygon_width test_polygon_width.cpp)
target_link_libraries(test_polygon_width PRIVATE geomlib)
add_test(NAME test_polygon_width COMMAND test_polygon_width)

add_executable(test_aabb_distance test_aabb_distance.cpp)
target_link_libraries(test_aabb_distance PRIVATE geomlib)
add_test(NAME test_aabb_distance COMMAND test_aabb_distance)

add_executable(test_ray_obb test_ray_obb.cpp)
target_link_libraries(test_ray_obb PRIVATE geomlib)
add_test(NAME test_ray_obb COMMAND test_ray_obb)

add_executable(test_tri_tri_distance test_tri_tri_distance.cpp)
target_link_libraries(test_tri_tri_distance PRIVATE geomlib)
add_test(NAME test_tri_tri_distance COMMAND test_tri_tri_distance)

add_executable(test_plane_sphere test_plane_sphere.cpp)
target_link_libraries(test_plane_sphere PRIVATE geomlib)
add_test(NAME test_plane_sphere COMMAND test_plane_sphere)

add_executable(test_circle_line test_circle_line.cpp)
target_link_libraries(test_circle_line PRIVATE geomlib)
add_test(NAME test_circle_line COMMAND test_circle_line)

add_executable(test_ray_tetra test_ray_tetra.cpp)
target_link_libraries(test_ray_tetra PRIVATE geomlib)
add_test(NAME test_ray_tetra COMMAND test_ray_tetra)

add_executable(test_polyline_closest test_polyline_closest.cpp)
target_link_libraries(test_polyline_closest PRIVATE geomlib)
add_test(NAME test_polyline_closest COMMAND test_polyline_closest)

add_executable(test_capsule_aabb test_capsule_aabb.cpp)
target_link_libraries(test_capsule_aabb PRIVATE geomlib)
add_test(NAME test_capsule_aabb COMMAND test_capsule_aabb)

add_executable(test_polygon_distance2d test_polygon_distance2d.cpp)
target_link_libraries(test_polygon_distance2d PRIVATE geomlib)
add_test(NAME test_polygon_distance2d COMMAND test_polygon_distance2d)

add_executable(test_tri_tetra test_tri_tetra.cpp)
target_link_libraries(test_tri_tetra PRIVATE geomlib)
add_test(NAME test_tri_tetra COMMAND test_tri_tetra)

add_executable(test_sphere_capsule test_sphere_capsule.cpp)
target_link_libraries(test_sphere_capsule PRIVATE geomlib)
add_test(NAME test_sphere_capsule COMMAND test_sphere_capsule)

add_executable(test_tetra_tetra test_tetra_tetra.cpp)
target_link_libraries(test_tetra_tetra PRIVATE geomlib)
add_test(NAME test_tetra_tetra COMMAND test_tetra_tetra)

add_executable(test_sphere_obb test_sphere_obb.cpp)
target_link_libraries(test_sphere_obb PRIVATE geomlib)
add_test(NAME test_sphere_obb COMMAND test_sphere_obb)

add_executable(test_circle_aabb test_circle_aabb.cpp)
target_link_libraries(test_circle_aabb PRIVATE geomlib)
add_test(NAME test_circle_aabb COMMAND test_circle_aabb)

add_executable(test_capsule_obb test_capsule_obb.cpp)
target_link_libraries(test_capsule_obb PRIVATE geomlib)
add_test(NAME test_capsule_obb COMMAND test_capsule_obb)

add_executable(test_ray_cone test_ray_cone.cpp)
target_link_libraries(test_ray_cone PRIVATE geomlib)
add_test(NAME test_ray_cone COMMAND test_ray_cone)

add_executable(test_tri_obb test_tri_obb.cpp)
target_link_libraries(test_tri_obb PRIVATE geomlib)
add_test(NAME test_tri_obb COMMAND test_tri_obb)

add_executable(test_aabb_obb test_aabb_obb.cpp)
target_link_libraries(test_aabb_obb PRIVATE geomlib)
add_test(NAME test_aabb_obb COMMAND test_aabb_obb)

add_executable(test_obb_obb test_obb_obb.cpp)
target_link_libraries(test_obb_obb PRIVATE geomlib)
add_test(NAME test_obb_obb COMMAND test_obb_obb)

add_executable(test_ellipse test_ellipse.cpp)
target_link_libraries(test_ellipse PRIVATE geomlib)
add_test(NAME test_ellipse COMMAND test_ellipse)

add_executable(test_capsule_triangle test_capsule_triangle.cpp)
target_link_libraries(test_capsule_triangle PRIVATE geomlib)
add_test(NAME test_capsule_triangle COMMAND test_capsule_triangle)

add_executable(test_ellipse_closest test_ellipse_closest.cpp)
target_link_libraries(test_ellipse_closest PRIVATE geomlib)
add_test(NAME test_ellipse_closest COMMAND test_ellipse_closest)

add_executable(test_ellipse_measures test_ellipse_measures.cpp)
target_link_libraries(test_ellipse_measures PRIVATE geomlib)
add_test(NAME test_ellipse_measures COMMAND test_ellipse_measures)

add_executable(test_circle_triangle test_circle_triangle.cpp)
target_link_libraries(test_circle_triangle PRIVATE geomlib)
add_test(NAME test_circle_triangle COMMAND test_circle_triangle)

add_executable(test_line2d test_line2d.cpp)
target_link_libraries(test_line2d PRIVATE geomlib)
add_test(NAME test_line2d COMMAND test_line2d)

add_executable(test_arc_length test_arc_length.cpp)
target_link_libraries(test_arc_length PRIVATE geomlib)
add_test(NAME test_arc_length COMMAND test_arc_length)

add_executable(test_circle_area test_circle_area.cpp)
target_link_libraries(test_circle_area PRIVATE geomlib)
add_test(NAME test_circle_area COMMAND test_circle_area)

add_executable(test_polygon_intersection_area test_polygon_intersection_area.cpp)
target_link_libraries(test_polygon_intersection_area PRIVATE geomlib)
add_test(NAME test_polygon_intersection_area COMMAND test_polygon_intersection_area)

add_executable(test_circle_ellipse test_circle_ellipse.cpp)
target_link_libraries(test_circle_ellipse PRIVATE geomlib)
add_test(NAME test_circle_ellipse COMMAND test_circle_ellipse)

add_executable(test_geodesic_mmp test_geodesic_mmp.cpp)
target_link_libraries(test_geodesic_mmp PRIVATE geomlib)
add_test(NAME test_geodesic_mmp COMMAND test_geodesic_mmp)

add_executable(test_metric test_metric.cpp)
target_link_libraries(test_metric PRIVATE geomlib)
add_test(NAME test_metric COMMAND test_metric)
