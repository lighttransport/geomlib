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

# Strict C11 for C sources (no GNU extensions). The web build mirrors these
# settings in web/CMakeLists.txt.
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(GEOMLIB_STRICT_C11 "Compile the C core with -pedantic -Werror" OFF)

# Shared source list, used by both the native and Emscripten builds.
set(GEOMLIB_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/geomlib_sources.cmake)

# Build options
option(GEOMLIB_BUILD_TESTS "Build tests" ON)
option(GEOMLIB_BUILD_BENCHMARKS "Build benchmark tools" OFF)
# GEOMLIB_BUILD_EXAMPLES was removed: no example targets exist.

# No RTTI, no exceptions
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti -fno-exceptions")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GR- /EHsc")
endif()

# Create library from the shared source list.
add_library(geomlib ${GEOMLIB_SOURCES})

if(GEOMLIB_STRICT_C11)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(geomlib PRIVATE -pedantic -Wall -Wextra -Werror)
    endif()
endif()

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

# Benchmarks
if(GEOMLIB_BUILD_BENCHMARKS)
    add_executable(benchmark_sdf scripts/benchmark_sdf.c)
    target_link_libraries(benchmark_sdf PRIVATE geomlib)
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
