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

# Create build directory
mkdir -p "$PROJECT_ROOT/$BUILD_DIR"
cd "$PROJECT_ROOT/$BUILD_DIR"

# Prepare CMake command
CMAKE_CMD="cmake"
CMAKE_CMD="$CMAKE_CMD -DCMAKE_BUILD_TYPE=$BUILD_TYPE"
CMAKE_CMD="$CMAKE_CMD -DGEOMLIB_BUILD_TESTS=$BUILD_TESTS"
CMAKE_CMD="$CMAKE_CMD -DGEOMLIB_BUILD_EXAMPLES=$BUILD_EXAMPLES"

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
echo "To build:         cd $BUILD_DIR && make"
echo "To test:          cd $BUILD_DIR && ctest"
echo "To install:       cd $BUILD_DIR && make install"
echo "==================================="
