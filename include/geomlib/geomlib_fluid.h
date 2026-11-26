/*
 * geomlib - Grid-based incompressible fluid solver
 *
 * Semi-Lagrangian advection followed by a pressure-projection step that
 * removes divergence from a velocity field on a regular 3D grid.
 */

#ifndef GEOMLIB_FLUID_H
#define GEOMLIB_FLUID_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t nx, ny, nz;
    float dx;
    float *u, *v, *w;      /* velocity components at cell centres */
    float *pressure;        /* last solved pressure field */
    float *divergence;      /* last divergence field */
} geomlib_fluid;

geomlib_fluid geomlib_fluid_create(uint32_t nx, uint32_t ny, uint32_t nz, float dx);
void geomlib_fluid_free(geomlib_fluid *f);

/* Remove divergence: solve Laplacian(p) = div(u)/dt and subtract grad(p). */
int geomlib_fluid_project(geomlib_fluid *f, float dt);

/* Semi-Lagrangian advection of the velocity field. */
int geomlib_fluid_advect(geomlib_fluid *f, float dt);

/* RMS divergence of the current field (useful for verification). */
float geomlib_fluid_rms_divergence(const geomlib_fluid *f);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_FLUID_H */
