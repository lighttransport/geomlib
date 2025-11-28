cmake_minimum_required(VERSION 3.14)

project(geomlib-wasm
    VERSION 0.1.0
    DESCRIPTION "geomlib WebAssembly build"
    LANGUAGES C CXX
)

# Check that we're using Emscripten
if(NOT EMSCRIPTEN)
    message(FATAL_ERROR "This CMakeLists.txt is for Emscripten build only. Use: emcmake cmake ..")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Source files from parent directory
set(GEOMLIB_SOURCES
    ../src/geomlib_linalg.c
    ../src/geomlib_solver.c
    ../src/geomlib_mesh.c
    ../src/geomlib_random.c
    ../src/geomlib_montecarlo.c
    ../src/geomlib_bvh.c
    ../src/geomlib_collision.c
    ../src/geomlib_surface.c
    ../src/geomlib_healing.c
    ../src/geomlib_remeshing.c
    ../src/geomlib_parameterization.c
    ../src/geomlib_curvature.c
    ../src/geomlib_geodesic.c
    ../src/geomlib_decimation.c
    ../src/geomlib_smoothing.c
    ../src/geomlib_quality.c
    ../src/geomlib_subdivision.c
    ../src/geomlib_atlas.c
    ../src/geomlib_tetramesh.c
    ../src/geomlib_cfd.c
    ../src/geomlib_optimize.c
    ../src/geomlib_topology.c
    ../src/geomlib_statistics.c
)

# WebAssembly binding
set(BINDING_SOURCE
    binding.cc
)

# Create wasm library
add_executable(geomlib-wasm ${GEOMLIB_SOURCES} ${BINDING_SOURCE})

target_include_directories(geomlib-wasm PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# Emscripten-specific settings
set(EMSCRIPTEN_LINK_FLAGS
    "-lembind"
    "--bind"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sMODULARIZE=1"
    "-sEXPORT_NAME=createGeomlibModule"
    "-sEXPORT_ES6=1"
    "-sENVIRONMENT=web"
    "-sMALLOC=emmalloc"
    "-sSTACK_SIZE=5MB"
    "-sINITIAL_MEMORY=16MB"
    "-sEXPORTED_RUNTIME_METHODS=['ccall','cwrap']"
    "-sNO_DISABLE_EXCEPTION_CATCHING"
)

# Add optimization flags for release
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    list(APPEND EMSCRIPTEN_LINK_FLAGS
        "-O3"
        "-flto"
    )
else()
    list(APPEND EMSCRIPTEN_LINK_FLAGS
        "-O0"
        "-g"
        "-sASSERTIONS=1"
    )
endif()

target_link_options(geomlib-wasm PRIVATE ${EMSCRIPTEN_LINK_FLAGS})

# Set output name
set_target_properties(geomlib-wasm PROPERTIES
    OUTPUT_NAME "geomlib"
    SUFFIX ".js"
)

# Copy output to demo directory after build
add_custom_command(TARGET geomlib-wasm POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/geomlib.js
        ${CMAKE_CURRENT_SOURCE_DIR}/demo/public/geomlib.js
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/geomlib.wasm
        ${CMAKE_CURRENT_SOURCE_DIR}/demo/public/geomlib.wasm
    COMMENT "Copying wasm files to demo/public/"
)
