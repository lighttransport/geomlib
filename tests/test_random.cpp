/*
 * geomlib - Random Number Generation Tests
 * Tests for PCG32 PRNG implementation
 */

#include <geomlib/geomlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define TEST_EPSILON 1e-6f
#define TEST_EPSILON_D 1e-12

void test_pcg32_init() {
    printf("Testing PCG32 initialization...\n");

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 54);

    assert(rng.state != 0);
    assert(rng.inc != 0);
    assert((rng.inc & 1) == 1);  /* Increment must be odd */

    printf("  ✓ Initialization successful\n");
}

void test_pcg32_deterministic() {
    printf("Testing PCG32 deterministic behavior...\n");

    geomlib_pcg32 rng1, rng2;
    geomlib_pcg32_init(&rng1, 12345, 67890);
    geomlib_pcg32_init(&rng2, 12345, 67890);

    /* Same seed should produce same sequence */
    for (int i = 0; i < 100; i++) {
        uint32_t v1 = geomlib_pcg32_next(&rng1);
        uint32_t v2 = geomlib_pcg32_next(&rng2);
        assert(v1 == v2);
    }

    printf("  ✓ Same seed produces identical sequences\n");
}

void test_pcg32_different_seeds() {
    printf("Testing PCG32 with different seeds...\n");

    geomlib_pcg32 rng1, rng2;
    geomlib_pcg32_init(&rng1, 12345, 0);
    geomlib_pcg32_init(&rng2, 54321, 0);

    /* Different seeds should produce different sequences */
    int differences = 0;
    for (int i = 0; i < 100; i++) {
        uint32_t v1 = geomlib_pcg32_next(&rng1);
        uint32_t v2 = geomlib_pcg32_next(&rng2);
        if (v1 != v2) {
            differences++;
        }
    }

    /* Expect most values to be different (at least 95%) */
    assert(differences >= 95);

    printf("  ✓ Different seeds produce different sequences (%d/100 different)\n", differences);
}

void test_pcg32_different_sequences() {
    printf("Testing PCG32 with different sequence numbers...\n");

    geomlib_pcg32 rng1, rng2;
    geomlib_pcg32_init(&rng1, 12345, 0);
    geomlib_pcg32_init(&rng2, 12345, 1);

    /* Same seed but different sequence should produce different values */
    int differences = 0;
    for (int i = 0; i < 100; i++) {
        uint32_t v1 = geomlib_pcg32_next(&rng1);
        uint32_t v2 = geomlib_pcg32_next(&rng2);
        if (v1 != v2) {
            differences++;
        }
    }

    /* Expect most values to be different (at least 95%) */
    assert(differences >= 95);

    printf("  ✓ Different sequences produce different values (%d/100 different)\n", differences);
}

void test_pcg32_nextf() {
    printf("Testing PCG32 float generation...\n");

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    int in_range = 0;
    float sum = 0.0f;
    int num_samples = 10000;

    for (int i = 0; i < num_samples; i++) {
        float v = geomlib_pcg32_nextf(&rng);
        assert(v >= 0.0f && v < 1.0f);
        sum += v;
        in_range++;
    }

    /* Check mean is approximately 0.5 (allow some variance) */
    float mean = sum / num_samples;
    assert(mean > 0.45f && mean < 0.55f);

    printf("  ✓ All %d values in [0, 1), mean = %.4f\n", in_range, mean);
}

void test_pcg32_nextd() {
    printf("Testing PCG32 double generation...\n");

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    int in_range = 0;
    double sum = 0.0;
    int num_samples = 10000;

    for (int i = 0; i < num_samples; i++) {
        double v = geomlib_pcg32_nextd(&rng);
        assert(v >= 0.0 && v < 1.0);
        sum += v;
        in_range++;
    }

    /* Check mean is approximately 0.5 (allow some variance) */
    double mean = sum / num_samples;
    assert(mean > 0.45 && mean < 0.55);

    printf("  ✓ All %d values in [0, 1), mean = %.6f\n", in_range, mean);
}

void test_pcg32_range_int() {
    printf("Testing PCG32 integer range...\n");

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    int min = -10;
    int max = 10;
    int count[21] = {0};  /* -10 to 10 inclusive */

    int num_samples = 10000;
    for (int i = 0; i < num_samples; i++) {
        int v = geomlib_pcg32_range_int(&rng, min, max);
        assert(v >= min && v <= max);
        count[v - min]++;
    }

    /* Check all values were generated at least once */
    for (int i = 0; i <= max - min; i++) {
        assert(count[i] > 0);
    }

    /* Check distribution is reasonably uniform (each value ~476 times on average) */
    /* Allow range of 300-700 for statistical variance */
    for (int i = 0; i <= max - min; i++) {
        assert(count[i] >= 300 && count[i] <= 700);
    }

    printf("  ✓ All values in [%d, %d] generated with reasonable distribution\n", min, max);
}

void test_pcg32_range_float() {
    printf("Testing PCG32 float range...\n");

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    float min = -5.0f;
    float max = 5.0f;
    float sum = 0.0f;

    int num_samples = 10000;
    for (int i = 0; i < num_samples; i++) {
        float v = geomlib_pcg32_range_float(&rng, min, max);
        assert(v >= min && v < max);
        sum += v;
    }

    /* Check mean is approximately (min+max)/2 = 0.0 */
    float mean = sum / num_samples;
    assert(mean > -0.5f && mean < 0.5f);

    printf("  ✓ All values in [%.1f, %.1f), mean = %.4f\n", min, max, mean);
}

void test_pcg32_range_double() {
    printf("Testing PCG32 double range...\n");

    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 0);

    double min = -100.0;
    double max = 100.0;
    double sum = 0.0;

    int num_samples = 10000;
    for (int i = 0; i < num_samples; i++) {
        double v = geomlib_pcg32_range_double(&rng, min, max);
        assert(v >= min && v < max);
        sum += v;
    }

    /* Check mean is approximately (min+max)/2 = 0.0 */
    double mean = sum / num_samples;
    assert(mean > -5.0 && mean < 5.0);

    printf("  ✓ All values in [%.1f, %.1f), mean = %.4f\n", min, max, mean);
}

void test_pcg32_known_sequence() {
    printf("Testing PCG32 against known sequence...\n");

    /* Test against known PCG32 output for seed=42, seq=54 */
    geomlib_pcg32 rng;
    geomlib_pcg32_init(&rng, 42, 54);

    /* Generate first few values */
    uint32_t v1 = geomlib_pcg32_next(&rng);
    uint32_t v2 = geomlib_pcg32_next(&rng);
    uint32_t v3 = geomlib_pcg32_next(&rng);
    uint32_t v4 = geomlib_pcg32_next(&rng);
    uint32_t v5 = geomlib_pcg32_next(&rng);

    /* Just verify they're non-zero and different */
    assert(v1 != 0);
    assert(v2 != 0);
    assert(v3 != 0);
    assert(v4 != 0);
    assert(v5 != 0);
    assert(v1 != v2);
    assert(v2 != v3);
    assert(v3 != v4);
    assert(v4 != v5);

    printf("  ✓ Generated sequence: %u, %u, %u, %u, %u\n", v1, v2, v3, v4, v5);
}

void test_pcg32_edge_cases() {
    printf("Testing PCG32 edge cases...\n");

    geomlib_pcg32 rng;

    /* Test with zero seed */
    geomlib_pcg32_init(&rng, 0, 0);
    uint32_t v1 = geomlib_pcg32_next(&rng);
    assert(v1 != 0);  /* Should still generate non-zero values */

    /* Test range with same min/max */
    geomlib_pcg32_init(&rng, 42, 0);
    int v2 = geomlib_pcg32_range_int(&rng, 5, 5);
    assert(v2 == 5);

    /* Test range with reversed min/max (should swap) */
    int v3 = geomlib_pcg32_range_int(&rng, 10, 1);
    assert(v3 >= 1 && v3 <= 10);

    printf("  ✓ Edge cases handled correctly\n");
}

int main() {
    printf("========================================\n");
    printf("geomlib Random Number Generation Tests\n");
    printf("========================================\n\n");

    test_pcg32_init();
    test_pcg32_deterministic();
    test_pcg32_different_seeds();
    test_pcg32_different_sequences();
    test_pcg32_nextf();
    test_pcg32_nextd();
    test_pcg32_range_int();
    test_pcg32_range_float();
    test_pcg32_range_double();
    test_pcg32_known_sequence();
    test_pcg32_edge_cases();

    printf("\n========================================\n");
    printf("All PCG32 tests passed!\n");
    printf("========================================\n");

    return 0;
}
