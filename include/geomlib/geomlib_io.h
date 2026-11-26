/*
 * geomlib - Mesh file I/O
 *
 * All entry points operate on in-memory buffers so they work identically in
 * native builds and inside WebAssembly (where there is no real filesystem).
 * File-path wrappers can be added by callers using their platform's FS API.
 *
 * Supported formats: OBJ, OFF, ASCII PLY, ASCII/binary STL.
 */

#ifndef GEOMLIB_IO_H
#define GEOMLIB_IO_H

#include "geomlib_mesh.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define GEOMLIB_IO_SUCCESS        0
#define GEOMLIB_IO_ERROR_MEMORY  (-1)
#define GEOMLIB_IO_ERROR_INVALID (-2)
#define GEOMLIB_IO_ERROR_UNSUPPORTED (-3)

/* Guess a format name ("obj", "off", "ply", "stl") from the buffer contents.
 * Returns NULL when unknown. */
const char *geomlib_mesh_detect_format(const void *data, size_t size);

/* Load a mesh from a memory buffer. format may be NULL for auto-detection.
 * On success the mesh is initialized and filled; on failure it is left empty. */
int geomlib_mesh_load_mem(geomlib_trimesh *mesh, const char *format,
                          const void *data, size_t size);

/* Serialize a mesh into a newly malloc'd buffer. The caller owns *out_data
 * and must free() it. */
int geomlib_mesh_save_mem(const geomlib_trimesh *mesh, const char *format,
                          void **out_data, size_t *out_size);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_IO_H */
