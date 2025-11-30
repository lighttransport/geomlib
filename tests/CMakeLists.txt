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
