/*
 * geomlib - 2D polygon rasterization
 *
 * Samples a simple polygon onto a regular grid, producing an inside mask and
 * optional signed distances at cell centres.
 */

#ifndef GEOMLIB_RASTER2D_H
#define GEOMLIB_RASTER2D_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t nx, ny;
    geomlib_vec2 origin;
    float cell_size;
    uint8_t *mask;     /* nx*ny, 1 when the cell centre is inside */
    float *distance;   /* optional signed distance, NULL when not requested */
} geomlib_raster2d;

geomlib_raster2d geomlib_raster2d_create(void);
void geomlib_raster2d_free(geomlib_raster2d *r);

/* Rasterize `polygon` into an nx*ny grid with margins. */
int geomlib_rasterize_polygon2d(const geomlib_vec2 *polygon, size_t n,
                                uint32_t nx, uint32_t ny, float margin,
                                int with_distance, geomlib_raster2d *out);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_RASTER2D_H */
