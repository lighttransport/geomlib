cmake_minimum_required(VERSION 3.14)

project(geomlib
    VERSION 0.1.0
    DESCRIPTION "A portable geometry processing library"
    LANGUAGES C CXX
)

# C++17 required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# C11 for C sources
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Build options
option(GEOMLIB_BUILD_TESTS "Build tests" ON)
option(GEOMLIB_BUILD_EXAMPLES "Build examples" OFF)

# No RTTI, no exceptions
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti -fno-exceptions")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GR- /EHsc")
endif()

# Library source files
set(GEOMLIB_SOURCES
    src/geomlib_linalg.c
    src/geomlib_solver.c
    src/geomlib_mesh.c
    src/geomlib_random.c
    src/geomlib_montecarlo.c
    src/geomlib_bvh.c
    src/geomlib_collision.c
    src/geomlib_surface.c
    src/geomlib_healing.c
    src/geomlib_remeshing.c
    src/geomlib_parameterization.c
)

# Create library
add_library(geomlib ${GEOMLIB_SOURCES})

target_include_directories(geomlib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Link math library on Unix
if(UNIX)
    target_link_libraries(geomlib PUBLIC m)
endif()

# Tests
if(GEOMLIB_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Installation
include(GNUInstallDirs)

install(TARGETS geomlib
    EXPORT geomlibTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/geomlib
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT geomlibTargets
    FILE geomlibTargets.cmake
    NAMESPACE geomlib::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/geomlib
)
