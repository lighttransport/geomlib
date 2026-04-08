/*
 * geomlib - A portable geometry processing library
 * MIT License - Copyright (c) 2025 Light Transport Entertainment, Inc.
 *
 * Custom allocator interface for freestanding and hosted environments.
 *
 * Design:
 * - Freestanding C11 compatible (only uses stddef.h, stdint.h)
 * - Global allocator set once at startup; all library allocations go through it
 * - Default allocator uses stdlib malloc/free in hosted environments
 * - Users must call geomlib_set_allocator() before any library function
 *   in freestanding mode (GEOMLIB_FREESTANDING defined)
 */

#ifndef GEOMLIB_ALLOCATOR_H
#define GEOMLIB_ALLOCATOR_H

/* Freestanding-safe includes only */
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------------------
 * Allocator callback interface
 * -------------------------------------------------------------------------- */

/**
 * Custom allocator structure.
 *
 * Users provide function pointers for memory allocation, reallocation,
 * and deallocation. An opaque user_data pointer is passed to every callback.
 *
 * Requirements:
 * - malloc_fn must return a pointer aligned to at least _Alignof(max_align_t)
 *   or NULL on failure. When size is 0, behavior is implementation-defined
 *   (may return NULL or a unique pointer).
 * - realloc_fn must behave like C realloc: if ptr is NULL, acts like malloc.
 *   Returns NULL on failure without freeing the original pointer.
 * - free_fn must accept NULL (no-op) and free memory allocated by malloc_fn
 *   or realloc_fn.
 */
typedef struct geomlib_allocator {
    void* (*malloc_fn)(size_t size, void *user_data);
    void* (*realloc_fn)(void *ptr, size_t new_size, void *user_data);
    void  (*free_fn)(void *ptr, void *user_data);
    void  *user_data;
} geomlib_allocator;

/* --------------------------------------------------------------------------
 * Global allocator management
 * -------------------------------------------------------------------------- */

/**
 * Set the global allocator used by all library functions.
 *
 * Must be called before any other geomlib function that allocates memory.
 * Pass NULL to reset to the default allocator (hosted mode only).
 *
 * Thread safety: Not thread-safe. Set once at program startup before
 * any concurrent library use.
 */
void geomlib_set_allocator(const geomlib_allocator *allocator);

/**
 * Get the current global allocator.
 *
 * Returns the allocator set by geomlib_set_allocator(), or the default
 * stdlib-based allocator if none was set (hosted mode only).
 * Returns NULL if no allocator is set in freestanding mode.
 */
const geomlib_allocator* geomlib_get_allocator(void);

#ifndef GEOMLIB_FREESTANDING
/**
 * Get the default allocator (uses stdlib malloc/realloc/free).
 * Only available in hosted environments.
 */
const geomlib_allocator* geomlib_default_allocator(void);
#endif

/* --------------------------------------------------------------------------
 * Allocation functions (used internally by the library)
 * -------------------------------------------------------------------------- */

/**
 * Allocate memory through the global allocator.
 * Returns NULL on failure or if no allocator is set.
 */
void* geomlib_malloc(size_t size);

/**
 * Allocate zero-initialized memory through the global allocator.
 * Returns NULL on failure, on overflow (count * size), or if no allocator is set.
 */
void* geomlib_calloc(size_t count, size_t size);

/**
 * Reallocate memory through the global allocator.
 * Returns NULL on failure without freeing the original pointer.
 */
void* geomlib_realloc(void *ptr, size_t new_size);

/**
 * Free memory through the global allocator.
 * Accepts NULL (no-op).
 */
void geomlib_free(void *ptr);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_ALLOCATOR_H */
