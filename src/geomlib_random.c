/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Random Number Generation - PCG32 Implementation
 */

#include "geomlib/geomlib_random.h"
#include <string.h>

/* PCG32 constants */
#define GEOMLIB_PCG32_MULT 0x5851f42d4c957f2dULL
#define GEOMLIB_PCG32_DEFAULT_INC 0x14057b7ef767814fULL

/* ============================================================================
 * Initialization
 * ============================================================================ */

void geomlib_pcg32_init(geomlib_pcg32 *rng, uint64_t seed, uint64_t seq) {
    /* Initialize state to zero */
    rng->state = 0U;

    /* Set increment (must be odd) */
    rng->inc = (seq << 1u) | 1u;

    /* Advance state once */
    geomlib_pcg32_next(rng);

    /* Mix in the seed */
    rng->state += seed;

    /* Advance state again */
    geomlib_pcg32_next(rng);
}

/* ============================================================================
 * Core Generator
 * ============================================================================ */

uint32_t geomlib_pcg32_next(geomlib_pcg32 *rng) {
    uint64_t oldstate = rng->state;

    /* Advance internal state (LCG) */
    rng->state = oldstate * GEOMLIB_PCG32_MULT + rng->inc;

    /* Calculate output function (XSH-RR):
     * XOR high and low bits, then shift right by some amount.
     * The shift amount comes from the highest bits of oldstate.
     * Finally, rotate right by the rotation amount.
     */
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = (uint32_t)(oldstate >> 59u);

    /* Rotate right: (x >> r) | (x << ((-r) & 31)) */
    return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
}

/* ============================================================================
 * Floating Point Generation
 * ============================================================================ */

float geomlib_pcg32_nextf(geomlib_pcg32 *rng) {
    /* Generate uniform float in [0, 1) using IEEE 754 bit manipulation.
     *
     * Method: Generate a float in [1, 2) by setting exponent bits to 127
     * (which represents 2^0 = 1 in biased exponent), then place random
     * bits in the mantissa. Subtract 1.0 to get [0, 1).
     *
     * This avoids expensive floating-point division.
     */
    union {
        uint32_t u;
        float f;
    } x;

    /* Take upper 23 bits for mantissa, combine with exponent for 1.0 */
    x.u = (geomlib_pcg32_next(rng) >> 9) | 0x3f800000u;

    /* Subtract 1.0 to get [0, 1) */
    return x.f - 1.0f;
}

double geomlib_pcg32_nextd(geomlib_pcg32 *rng) {
    /* Generate uniform double in [0, 1) using IEEE 754 bit manipulation.
     *
     * Method: Generate a double in [1, 2) by setting exponent bits to 1023
     * (which represents 2^0 = 1 in biased exponent), then place random
     * bits in the mantissa. Subtract 1.0 to get [0, 1).
     *
     * This requires 52 bits of mantissa, so we use a single 32-bit random
     * number shifted to the upper bits (sacrificing some lower precision
     * for speed).
     */
    union {
        uint64_t u;
        double d;
    } x;

    /* Take 32 bits, shift to upper mantissa bits, combine with exponent for 1.0 */
    x.u = ((uint64_t)geomlib_pcg32_next(rng) << 20) | 0x3ff0000000000000ULL;

    /* Subtract 1.0 to get [0, 1) */
    return x.d - 1.0;
}

/* ============================================================================
 * Range Generation
 * ============================================================================ */

int32_t geomlib_pcg32_range_int(geomlib_pcg32 *rng, int32_t min, int32_t max) {
    if (min > max) {
        int32_t temp = min;
        min = max;
        max = temp;
    }

    if (min == max) {
        return min;
    }

    /* Calculate range size (as unsigned to handle full int32 range) */
    uint32_t range = (uint32_t)(max - min) + 1u;

    /* Use rejection sampling to avoid modulo bias
     * This ensures uniform distribution across the range
     */
    uint32_t threshold = (~range + 1u) % range;
    uint32_t r;

    do {
        r = geomlib_pcg32_next(rng);
    } while (r < threshold);

    return min + (int32_t)(r % range);
}

float geomlib_pcg32_range_float(geomlib_pcg32 *rng, float min, float max) {
    return min + (max - min) * geomlib_pcg32_nextf(rng);
}

double geomlib_pcg32_range_double(geomlib_pcg32 *rng, double min, double max) {
    return min + (max - min) * geomlib_pcg32_nextd(rng);
}
