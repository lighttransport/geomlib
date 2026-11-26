#!/bin/bash

# CMake configuration script for geomlib

set -e

# Default values
BUILD_TYPE="Release"
BUILD_DIR="build"
BUILD_TESTS="ON"
BUILD_EXAMPLES="OFF"
INSTALL_PREFIX=""
EXTRA_CMAKE_ARGS=""
WASM_BUILD="OFF"
STRICT_C11="OFF"

# Print usage
usage() {
    cat << EOF
Usage: $0 [options]

Options:
    -h, --help              Show this help message
    -d, --debug             Build in Debug mode (default: Release)
    -b, --build-dir DIR     Build directory (default: build)
    -t, --no-tests          Disable building tests (default: enabled)
    -e, --examples          Enable building examples (default: disabled)
    -p, --prefix PATH       Installation prefix
    --wasm                  Configure the Emscripten WebAssembly build
    --strict-c11            Enforce strict C11 (-pedantic -Werror) in CI
    --cmake-args ARGS       Additional CMake arguments

Examples:
    $0                                  # Configure with defaults
    $0 --debug                          # Debug build
    $0 --build-dir build-debug --debug  # Debug build in custom directory
    $0 --prefix /usr/local              # Set install prefix
    $0 --examples                       # Enable examples

EOF
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            ;;
        -d|--debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        -b|--build-dir)
            BUILD_DIR="$2"
            shift 2
            ;;
        -t|--no-tests)
            BUILD_TESTS="OFF"
            shift
            ;;
        -e|--examples)
            BUILD_EXAMPLES="ON"
            shift
            ;;
        -p|--prefix)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        --wasm)
            WASM_BUILD="ON"
            shift
            ;;
        --strict-c11)
            STRICT_C11="ON"
            shift
            ;;
        --cmake-args)
            EXTRA_CMAKE_ARGS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "==================================="
echo "Configuring geomlib"
echo "==================================="
echo "Build type:       $BUILD_TYPE"
echo "Build directory:  $BUILD_DIR"
echo "Build tests:      $BUILD_TESTS"
echo "Build examples:   $BUILD_EXAMPLES"
if [ -n "$INSTALL_PREFIX" ]; then
    echo "Install prefix:   $INSTALL_PREFIX"
fi
echo "==================================="

# The Emscripten project lives in web/, so its build directory is web/build.
if [ "$WASM_BUILD" = "ON" ] && [ "$BUILD_DIR" = "build" ]; then
    BUILD_DIR="web/build"
fi

# Create build directory
mkdir -p "$PROJECT_ROOT/$BUILD_DIR"
cd "$PROJECT_ROOT/$BUILD_DIR"

# Prepare CMake command
if [ "$WASM_BUILD" = "ON" ]; then
    CMAKE_CMD="emcmake cmake"
else
    CMAKE_CMD="cmake"
fi
CMAKE_CMD="$CMAKE_CMD -DCMAKE_BUILD_TYPE=$BUILD_TYPE"
CMAKE_CMD="$CMAKE_CMD -DGEOMLIB_BUILD_TESTS=$BUILD_TESTS"
if [ "$BUILD_EXAMPLES" = "ON" ]; then
    echo "Warning: --examples is deprecated; no example targets exist." >&2
fi
CMAKE_CMD="$CMAKE_CMD -DGEOMLIB_STRICT_C11=$STRICT_C11"

if [ -n "$INSTALL_PREFIX" ]; then
    CMAKE_CMD="$CMAKE_CMD -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX"
fi

if [ -n "$EXTRA_CMAKE_ARGS" ]; then
    CMAKE_CMD="$CMAKE_CMD $EXTRA_CMAKE_ARGS"
fi

CMAKE_CMD="$CMAKE_CMD .."

# Run CMake
echo ""
echo "Running: $CMAKE_CMD"
echo ""

eval $CMAKE_CMD

echo ""
echo "==================================="
echo "Configuration complete!"
echo "==================================="
if [ "$WASM_BUILD" = "ON" ]; then
    echo "To build:         cd $BUILD_DIR && emmake make"
    echo "To test:          cd web && node smoke.mjs"
else
    echo "To build:         cd $BUILD_DIR && make"
    echo "To test:          cd $BUILD_DIR && ctest"
    echo "To install:       cd $BUILD_DIR && make install"
fi
echo "==================================="
