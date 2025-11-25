add_executable(test_linalg test_linalg.cpp)
target_link_libraries(test_linalg PRIVATE geomlib)
add_test(NAME test_linalg COMMAND test_linalg)

add_executable(test_solver test_solver.cpp)
target_link_libraries(test_solver PRIVATE geomlib)
add_test(NAME test_solver COMMAND test_solver)

add_executable(test_mesh test_mesh.cpp)
target_link_libraries(test_mesh PRIVATE geomlib)
add_test(NAME test_mesh COMMAND test_mesh)
