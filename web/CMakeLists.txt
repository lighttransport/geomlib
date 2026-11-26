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

# Strict C11 for the core, matching the native build.
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(GEOMLIB_STRICT_C11 "Compile the C core with -pedantic -Werror" OFF)

# Shared source list (identical to the native build).
set(GEOMLIB_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
include(${GEOMLIB_ROOT}/cmake/geomlib_sources.cmake)

# WebAssembly binding
set(BINDING_SOURCE
    binding.cc
)

# Create wasm library
add_executable(geomlib-wasm ${GEOMLIB_SOURCES} ${BINDING_SOURCE})

target_include_directories(geomlib-wasm PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

if(GEOMLIB_STRICT_C11)
    target_compile_options(geomlib-wasm PRIVATE -pedantic -Wall -Wextra -Werror)
endif()

# Emscripten-specific settings
set(EMSCRIPTEN_LINK_FLAGS
    "-lembind"
    "--bind"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sMODULARIZE=1"
    "-sEXPORT_NAME=createGeomlibModule"
    "-sEXPORT_ES6=1"
    "-sENVIRONMENT=web,node"
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
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/demo/public
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/geomlib.js
        ${CMAKE_CURRENT_SOURCE_DIR}/demo/public/geomlib.js
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/geomlib.wasm
        ${CMAKE_CURRENT_SOURCE_DIR}/demo/public/geomlib.wasm
    COMMENT "Copying wasm files to demo/public/"
)
