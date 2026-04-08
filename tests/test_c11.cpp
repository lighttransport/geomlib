/*
 * Tests for C11 allocator, RAII cleanup, safety utilities, and freestanding support.
 */

#include <cstdio>
#include <cstring>
#include <cstdlib>

extern "C" {
#include "geomlib/geomlib_allocator.h"
#include "geomlib/geomlib_c11.h"
#include "geomlib/geomlib_mesh.h"
#include "geomlib/geomlib_bvh.h"
#include "geomlib/geomlib_montecarlo.h"
}

static int test_count = 0;
static int pass_count = 0;

static void test(bool condition, const char* name) {
    test_count++;
    if (condition) {
        pass_count++;
        printf("[PASS] %s\n", name);
    } else {
        printf("[FAIL] %s\n", name);
    }
}

/* ======================================================================
 * Custom allocator for testing
 * ====================================================================== */

struct tracking_allocator_data {
    size_t alloc_count;
    size_t free_count;
    size_t total_allocated;
    size_t current_allocated;
};

static void* tracking_malloc(size_t size, void *user_data) {
    struct tracking_allocator_data *data = (struct tracking_allocator_data *)user_data;
    data->alloc_count++;
    data->total_allocated += size;
    data->current_allocated += size;
    /* Store size prefix for tracking frees */
    void *ptr = malloc(size + sizeof(size_t));
    if (!ptr) return NULL;
    *(size_t*)ptr = size;
    return (char*)ptr + sizeof(size_t);
}

static void* tracking_realloc(void *ptr, size_t new_size, void *user_data) {
    struct tracking_allocator_data *data = (struct tracking_allocator_data *)user_data;
    if (ptr == NULL) {
        return tracking_malloc(new_size, user_data);
    }
    void *real_ptr = (char*)ptr - sizeof(size_t);
    size_t old_size = *(size_t*)real_ptr;
    data->current_allocated -= old_size;
    data->current_allocated += new_size;
    data->total_allocated += new_size;
    void *new_ptr = realloc(real_ptr, new_size + sizeof(size_t));
    if (!new_ptr) {
        data->current_allocated += old_size;
        data->current_allocated -= new_size;
        data->total_allocated -= new_size;
        return NULL;
    }
    *(size_t*)new_ptr = new_size;
    return (char*)new_ptr + sizeof(size_t);
}

static void tracking_free(void *ptr, void *user_data) {
    struct tracking_allocator_data *data = (struct tracking_allocator_data *)user_data;
    if (ptr == NULL) return;
    data->free_count++;
    void *real_ptr = (char*)ptr - sizeof(size_t);
    size_t size = *(size_t*)real_ptr;
    data->current_allocated -= size;
    free(real_ptr);
}

/* ======================================================================
 * Error handler for testing
 * ====================================================================== */

struct error_log {
    int count;
    geomlib_error_level last_level;
    char last_message[256];
};

static void test_error_handler(geomlib_error_level level, const char *message,
                                const char *file, int line, void *user_data) {
    struct error_log *log = (struct error_log *)user_data;
    log->count++;
    log->last_level = level;
    if (message) {
        strncpy(log->last_message, message, sizeof(log->last_message) - 1);
        log->last_message[sizeof(log->last_message) - 1] = '\0';
    }
    (void)file;
    (void)line;
}

/* ======================================================================
 * Test: Default allocator
 * ====================================================================== */

static void test_default_allocator(void) {
    /* Reset to default */
    geomlib_set_allocator(NULL);

    const geomlib_allocator *alloc = geomlib_get_allocator();
    test(alloc != NULL, "default allocator is non-NULL");
    test(alloc == geomlib_default_allocator(), "default allocator matches");

    /* Test basic allocation */
    void *ptr = geomlib_malloc(128);
    test(ptr != NULL, "geomlib_malloc returns non-NULL");
    memset(ptr, 0xAB, 128); /* Should not crash */
    geomlib_free(ptr);

    /* Test calloc */
    int *arr = (int*)geomlib_calloc(10, sizeof(int));
    test(arr != NULL, "geomlib_calloc returns non-NULL");
    bool all_zero = true;
    for (int i = 0; i < 10; i++) {
        if (arr[i] != 0) all_zero = false;
    }
    test(all_zero, "geomlib_calloc zero-initializes memory");
    geomlib_free(arr);

    /* Test realloc */
    ptr = geomlib_malloc(64);
    test(ptr != NULL, "geomlib_malloc for realloc test");
    memset(ptr, 0xCD, 64);
    void *new_ptr = geomlib_realloc(ptr, 256);
    test(new_ptr != NULL, "geomlib_realloc returns non-NULL");
    /* Check first 64 bytes are preserved */
    unsigned char *bytes = (unsigned char*)new_ptr;
    bool preserved = true;
    for (int i = 0; i < 64; i++) {
        if (bytes[i] != 0xCD) preserved = false;
    }
    test(preserved, "geomlib_realloc preserves existing data");
    geomlib_free(new_ptr);

    /* Test free NULL (should not crash) */
    geomlib_free(NULL);
    test(true, "geomlib_free(NULL) does not crash");
}

/* ======================================================================
 * Test: Custom allocator
 * ====================================================================== */

static void test_custom_allocator(void) {
    struct tracking_allocator_data data = {0, 0, 0, 0};
    geomlib_allocator custom = {
        tracking_malloc,
        tracking_realloc,
        tracking_free,
        &data
    };
    geomlib_set_allocator(&custom);

    test(geomlib_get_allocator() == &custom, "custom allocator is set");

    /* Allocate and free */
    void *p1 = geomlib_malloc(100);
    test(p1 != NULL, "custom malloc returns non-NULL");
    test(data.alloc_count == 1, "custom alloc_count is 1");
    test(data.current_allocated == 100, "custom tracks 100 bytes allocated");

    void *p2 = geomlib_malloc(200);
    test(data.alloc_count == 2, "custom alloc_count is 2");
    test(data.current_allocated == 300, "custom tracks 300 bytes total");

    geomlib_free(p1);
    test(data.free_count == 1, "custom free_count is 1");
    test(data.current_allocated == 200, "custom tracks 200 bytes after free");

    geomlib_free(p2);
    test(data.free_count == 2, "custom free_count is 2");
    test(data.current_allocated == 0, "custom tracks 0 bytes after all freed");

    /* Reset to default */
    geomlib_set_allocator(NULL);
}

/* ======================================================================
 * Test: Custom allocator with library types
 * ====================================================================== */

static void test_custom_allocator_with_mesh(void) {
    struct tracking_allocator_data data = {0, 0, 0, 0};
    geomlib_allocator custom = {
        tracking_malloc,
        tracking_realloc,
        tracking_free,
        &data
    };
    geomlib_set_allocator(&custom);

    /* Create and use a trimesh through the custom allocator */
    geomlib_trimesh mesh;
    int ret = geomlib_trimesh_init(&mesh, 10, 10);
    test(ret == 0, "trimesh_init succeeds with custom allocator");
    test(data.alloc_count > 0, "custom allocator was called for trimesh_init");

    size_t allocs_after_init = data.alloc_count;

    /* Add some vertices */
    geomlib_vec3 v = {1.0f, 2.0f, 3.0f};
    geomlib_trimesh_add_vertex(&mesh, v, NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, v, NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, v, NULL, NULL);
    test(mesh.vertex_count == 3, "trimesh has 3 vertices");

    /* Free the mesh */
    geomlib_trimesh_free(&mesh);
    test(data.current_allocated == 0, "all mesh memory freed through custom allocator");
    test(data.free_count > 0, "custom free was called");

    /* Reset */
    geomlib_set_allocator(NULL);
}

/* ======================================================================
 * Test: RAII-like cleanup
 * ====================================================================== */

#if GEOMLIB_HAS_CLEANUP
static void test_raii_trimesh(void) {
    struct tracking_allocator_data data = {0, 0, 0, 0};
    geomlib_allocator custom = {
        tracking_malloc,
        tracking_realloc,
        tracking_free,
        &data
    };
    geomlib_set_allocator(&custom);

    {
        GEOMLIB_AUTO_TRIMESH geomlib_trimesh mesh;
        geomlib_trimesh_init(&mesh, 8, 8);

        geomlib_vec3 v = {1.0f, 0.0f, 0.0f};
        geomlib_trimesh_add_vertex(&mesh, v, NULL, NULL);
        test(data.alloc_count > 0, "RAII: allocations made for trimesh");
    }
    /* mesh should be automatically freed when it goes out of scope */
    test(data.current_allocated == 0, "RAII: trimesh automatically freed on scope exit");

    geomlib_set_allocator(NULL);
}

static void test_raii_polygon2d(void) {
    struct tracking_allocator_data data = {0, 0, 0, 0};
    geomlib_allocator custom = {
        tracking_malloc,
        tracking_realloc,
        tracking_free,
        &data
    };
    geomlib_set_allocator(&custom);

    {
        GEOMLIB_AUTO_POLYGON2D geomlib_polygon2d poly;
        geomlib_polygon2d_init(&poly, 4);
        geomlib_vec2 v = {1.0f, 0.0f};
        geomlib_polygon2d_add_vertex(&poly, v);
        test(data.alloc_count > 0, "RAII: allocations made for polygon2d");
    }
    test(data.current_allocated == 0, "RAII: polygon2d automatically freed on scope exit");

    geomlib_set_allocator(NULL);
}

static void test_raii_auto_free_ptr(void) {
    struct tracking_allocator_data data = {0, 0, 0, 0};
    geomlib_allocator custom = {
        tracking_malloc,
        tracking_realloc,
        tracking_free,
        &data
    };
    geomlib_set_allocator(&custom);

    {
        GEOMLIB_AUTO_FREE void *ptr = geomlib_malloc(512);
        test(ptr != NULL, "RAII: auto_free pointer allocated");
        test(data.current_allocated == 512, "RAII: 512 bytes tracked");
    }
    test(data.current_allocated == 0, "RAII: auto_free pointer freed on scope exit");

    geomlib_set_allocator(NULL);
}
#endif /* GEOMLIB_HAS_CLEANUP */

/* ======================================================================
 * Test: Safe arithmetic
 * ====================================================================== */

static void test_safe_arithmetic(void) {
    size_t result;

    /* Normal multiply */
    test(geomlib_safe_mul_size(100, 200, &result) == 0, "safe_mul: 100*200 succeeds");
    test(result == 20000, "safe_mul: 100*200 = 20000");

    /* Multiply by zero */
    test(geomlib_safe_mul_size(0, 100, &result) == 0, "safe_mul: 0*100 succeeds");
    test(result == 0, "safe_mul: 0*100 = 0");

    /* Overflow */
    test(geomlib_safe_mul_size((size_t)-1, 2, &result) != 0, "safe_mul: SIZE_MAX*2 overflows");

    /* Large multiply that overflows */
    size_t half_max = ((size_t)-1) / 2 + 1;
    test(geomlib_safe_mul_size(half_max, 3, &result) != 0, "safe_mul: large values overflow");

    /* Normal add */
    test(geomlib_safe_add_size(100, 200, &result) == 0, "safe_add: 100+200 succeeds");
    test(result == 300, "safe_add: 100+200 = 300");

    /* Add overflow */
    test(geomlib_safe_add_size((size_t)-1, 1, &result) != 0, "safe_add: SIZE_MAX+1 overflows");

    /* NULL result pointer */
    test(geomlib_safe_mul_size(1, 1, NULL) != 0, "safe_mul: NULL result returns error");
    test(geomlib_safe_add_size(1, 1, NULL) != 0, "safe_add: NULL result returns error");
}

/* ======================================================================
 * Test: Calloc overflow protection
 * ====================================================================== */

static void test_calloc_overflow(void) {
    /* This should fail due to overflow */
    void *ptr = geomlib_calloc((size_t)-1, (size_t)-1);
    test(ptr == NULL, "geomlib_calloc detects overflow and returns NULL");
}

/* ======================================================================
 * Test: Error handling
 * ====================================================================== */

static void test_error_handling(void) {
    struct error_log log;
    memset(&log, 0, sizeof(log));

    geomlib_set_error_handler(test_error_handler, &log);

    /* Report an error */
    GEOMLIB_ERR("test error message");
    test(log.count == 1, "error handler called once");
    test(log.last_level == GEOMLIB_ERROR_ERROR, "error level is ERROR");
    test(strcmp(log.last_message, "test error message") == 0, "error message matches");

    /* Report a warning */
    GEOMLIB_WARN("test warning");
    test(log.count == 2, "error handler called again for warning");
    test(log.last_level == GEOMLIB_ERROR_WARNING, "error level is WARNING");

    /* Disable error handler */
    geomlib_set_error_handler(NULL, NULL);
    GEOMLIB_ERR("should not be seen");
    test(log.count == 2, "no handler: error count unchanged");
}

/* ======================================================================
 * Test: Memory utilities (geomlib_memset, etc.)
 * ====================================================================== */

static void test_memory_utilities(void) {
    /* geomlib_memset */
    char buf[64];
    geomlib_memset(buf, 0xAA, sizeof(buf));
    bool all_aa = true;
    for (size_t i = 0; i < sizeof(buf); i++) {
        if ((unsigned char)buf[i] != 0xAA) all_aa = false;
    }
    test(all_aa, "geomlib_memset fills correctly");

    /* geomlib_memcpy */
    char src[32], dst[32];
    for (int i = 0; i < 32; i++) src[i] = (char)i;
    geomlib_memcpy(dst, src, 32);
    test(geomlib_memcmp(dst, src, 32) == 0, "geomlib_memcpy copies correctly");

    /* geomlib_memcmp */
    test(geomlib_memcmp("abc", "abc", 3) == 0, "geomlib_memcmp: equal strings");
    test(geomlib_memcmp("abc", "abd", 3) < 0, "geomlib_memcmp: less than");
    test(geomlib_memcmp("abd", "abc", 3) > 0, "geomlib_memcmp: greater than");

    /* geomlib_memmove (overlapping) */
    char overlap[16] = "Hello, World!";
    geomlib_memmove(overlap + 2, overlap, 5);
    test(overlap[2] == 'H' && overlap[3] == 'e' && overlap[4] == 'l', 
         "geomlib_memmove handles overlapping regions");
}

/* ======================================================================
 * Test: Allocator with full pipeline (BVH)
 * ====================================================================== */

static void test_allocator_bvh_pipeline(void) {
    struct tracking_allocator_data data = {0, 0, 0, 0};
    geomlib_allocator custom = {
        tracking_malloc,
        tracking_realloc,
        tracking_free,
        &data
    };
    geomlib_set_allocator(&custom);

    /* Create a simple mesh */
    geomlib_trimesh mesh;
    geomlib_trimesh_init(&mesh, 4, 2);

    geomlib_vec3 v0 = {0, 0, 0}, v1 = {1, 0, 0}, v2 = {0, 1, 0}, v3 = {1, 1, 0};
    geomlib_trimesh_add_vertex(&mesh, v0, NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, v1, NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, v2, NULL, NULL);
    geomlib_trimesh_add_vertex(&mesh, v3, NULL, NULL);
    geomlib_trimesh_add_triangle(&mesh, 0, 1, 2);
    geomlib_trimesh_add_triangle(&mesh, 1, 3, 2);

    /* Build BVH */
    geomlib_bvh bvh;
    geomlib_bvh_build(&bvh, &mesh, GEOMLIB_BVH_QUALITY_MEDIUM);

    size_t peak = data.current_allocated;
    test(peak > 0, "BVH pipeline: memory allocated through custom allocator");

    /* Trace a ray */
    geomlib_ray ray = geomlib_ray_create(
        (geomlib_vec3){0.25f, 0.25f, 1.0f},
        (geomlib_vec3){0.0f, 0.0f, -1.0f}
    );
    geomlib_ray_hit hit;
    bool did_hit = geomlib_bvh_trace_ray(&bvh, &ray, &hit);
    test(did_hit, "BVH pipeline: ray hits mesh through custom allocator");

    /* Cleanup */
    geomlib_bvh_free(&bvh);
    geomlib_trimesh_free(&mesh);
    test(data.current_allocated == 0, "BVH pipeline: all memory freed");

    geomlib_set_allocator(NULL);
}

/* ======================================================================
 * Main
 * ====================================================================== */

int main(void) {
    printf("=== C11 Allocator, RAII, and Safety Tests ===\n\n");

    printf("--- Default Allocator ---\n");
    test_default_allocator();

    printf("\n--- Custom Allocator ---\n");
    test_custom_allocator();

    printf("\n--- Custom Allocator with Mesh ---\n");
    test_custom_allocator_with_mesh();

#if GEOMLIB_HAS_CLEANUP
    printf("\n--- RAII Cleanup (trimesh) ---\n");
    test_raii_trimesh();

    printf("\n--- RAII Cleanup (polygon2d) ---\n");
    test_raii_polygon2d();

    printf("\n--- RAII Cleanup (auto_free ptr) ---\n");
    test_raii_auto_free_ptr();
#else
    printf("\n--- RAII Cleanup: SKIPPED (compiler does not support __attribute__((cleanup))) ---\n");
#endif

    printf("\n--- Safe Arithmetic ---\n");
    test_safe_arithmetic();

    printf("\n--- Calloc Overflow Protection ---\n");
    test_calloc_overflow();

    printf("\n--- Error Handling ---\n");
    test_error_handling();

    printf("\n--- Memory Utilities ---\n");
    test_memory_utilities();

    printf("\n--- Allocator with BVH Pipeline ---\n");
    test_allocator_bvh_pipeline();

    printf("\n=== Results: %d/%d tests passed ===\n", pass_count, test_count);

    return (pass_count == test_count) ? 0 : 1;
}
