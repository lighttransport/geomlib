/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Random Number Generation - PCG32
 * Core functions implemented in C for maximum portability
 *
 * Based on PCG family by Melissa O'Neill: https://www.pcg-random.org/
 * Implementation technique from: https://github.com/wjakob/pcg32
 */

#ifndef GEOMLIB_RANDOM_H
#define GEOMLIB_RANDOM_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * PCG32 Random Number Generator State
 * ============================================================================ */

/**
 * PCG32 generator state
 *
 * PCG (Permuted Congruential Generator) is a family of simple, fast,
 * space-efficient statistically good algorithms with minimal overhead.
 */
typedef struct geomlib_pcg32 {
    uint64_t state;  /* RNG state */
    uint64_t inc;    /* Stream selection constant (must be odd) */
} geomlib_pcg32;

/* ============================================================================
 * Initialization
 * ============================================================================ */

/**
 * Initialize PCG32 generator with seed and sequence
 *
 * @param rng       Pointer to PCG32 state
 * @param seed      Initial seed value
 * @param seq       Sequence selection (determines the stream)
 *
 * Different sequences produce independent random streams.
 */
void geomlib_pcg32_init(geomlib_pcg32 *rng, uint64_t seed, uint64_t seq);

/* ============================================================================
 * Random Number Generation
 * ============================================================================ */

/**
 * Generate random 32-bit unsigned integer
 *
 * @param rng       Pointer to PCG32 state
 * @return          Random 32-bit unsigned integer
 */
uint32_t geomlib_pcg32_next(geomlib_pcg32 *rng);

/**
 * Generate random float in [0, 1)
 *
 * Uses IEEE 754 bit manipulation for efficient conversion without division.
 *
 * @param rng       Pointer to PCG32 state
 * @return          Random float in [0, 1)
 */
float geomlib_pcg32_nextf(geomlib_pcg32 *rng);

/**
 * Generate random double in [0, 1)
 *
 * Uses IEEE 754 bit manipulation for efficient conversion without division.
 *
 * @param rng       Pointer to PCG32 state
 * @return          Random double in [0, 1)
 */
double geomlib_pcg32_nextd(geomlib_pcg32 *rng);

/**
 * Generate random 32-bit signed integer in range [min, max]
 *
 * @param rng       Pointer to PCG32 state
 * @param min       Minimum value (inclusive)
 * @param max       Maximum value (inclusive)
 * @return          Random integer in [min, max]
 */
int32_t geomlib_pcg32_range_int(geomlib_pcg32 *rng, int32_t min, int32_t max);

/**
 * Generate random float in range [min, max)
 *
 * @param rng       Pointer to PCG32 state
 * @param min       Minimum value (inclusive)
 * @param max       Maximum value (exclusive)
 * @return          Random float in [min, max)
 */
float geomlib_pcg32_range_float(geomlib_pcg32 *rng, float min, float max);

/**
 * Generate random double in range [min, max)
 *
 * @param rng       Pointer to PCG32 state
 * @param min       Minimum value (inclusive)
 * @param max       Maximum value (exclusive)
 * @return          Random double in [min, max)
 */
double geomlib_pcg32_range_double(geomlib_pcg32 *rng, double min, double max);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_RANDOM_H */
