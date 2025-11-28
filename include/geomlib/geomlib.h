/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Main header file - includes all geomlib components
 *
 * Features:
 * - Linear algebra (vectors, matrices, transformations)
 * - Linear system solvers (Gaussian elimination, LU, iterative methods)
 * - Mesh/Polygon data structures and operations
 *
 * Design principles:
 * - C++17 compatible
 * - Core functions implemented in C for maximum portability
 * - No RTTI required
 * - No C++ exceptions
 * - No external dependencies (no Eigen, etc.)
 */

#ifndef GEOMLIB_H
#define GEOMLIB_H

/* Version information */
#define GEOMLIB_VERSION_MAJOR 0
#define GEOMLIB_VERSION_MINOR 1
#define GEOMLIB_VERSION_PATCH 0

/* Core components */
#include "geomlib_linalg.h"
#include "geomlib_solver.h"
#include "geomlib_mesh.h"
#include "geomlib_random.h"
#include "geomlib_montecarlo.h"
#include "geomlib_bvh.h"
#include "geomlib_collision.h"
#include "geomlib_surface.h"
#include "geomlib_healing.h"
#include "geomlib_remeshing.h"
#include "geomlib_parameterization.h"
#include "geomlib_curvature.h"
#include "geomlib_geodesic.h"
#include "geomlib_decimation.h"
#include "geomlib_smoothing.h"
#include "geomlib_quality.h"
#include "geomlib_subdivision.h"
#include "geomlib_atlas.h"
#include "geomlib_tetramesh.h"
#include "geomlib_cfd.h"
#include "geomlib_optimize.h"
#include "geomlib_topology.h"
#include "geomlib_statistics.h"
#include "geomlib_holefill.h"
#include "geomlib_feature.h"

#endif /* GEOMLIB_H */
