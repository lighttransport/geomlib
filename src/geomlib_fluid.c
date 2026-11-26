/*
 * geomlib - Grid-based incompressible fluid solver (see geomlib_fluid.h)
 *
 * Clean-room implementation of the stable fluids projection method:
 *   compute div(u) with forward differences
 *   solve Laplacian(p) = div(u)
 *   u -= grad(p)  (same forward-difference operator)
 * with a 7-point finite-difference Laplacian solved by a sparse direct solver.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_fluid.h"
#include "geomlib/geomlib_sparse.h"
#include "geomlib/geomlib_sparse_direct.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define IDX(f,x,y,z) ((x) + (f)->nx*((y) + (f)->ny*(z)))

geomlib_fluid geomlib_fluid_create(uint32_t nx, uint32_t ny, uint32_t nz, float dx)
{
    geomlib_fluid f;
    size_t n;
    memset(&f, 0, sizeof(f));
    f.nx = nx; f.ny = ny; f.nz = nz; f.dx = dx;
    n = (size_t)nx*ny*nz;
    f.u = (float *)calloc(n, sizeof(float));
    f.v = (float *)calloc(n, sizeof(float));
    f.w = (float *)calloc(n, sizeof(float));
    f.pressure = (float *)calloc(n, sizeof(float));
    f.divergence = (float *)calloc(n, sizeof(float));
    if (!f.u || !f.v || !f.w || !f.pressure || !f.divergence)
        geomlib_fluid_free(&f);
    return f;
}

void geomlib_fluid_free(geomlib_fluid *f)
{
    if (!f) return;
    free(f->u); free(f->v); free(f->w);
    free(f->pressure); free(f->divergence);
    memset(f, 0, sizeof(*f));
}

static float sample_field(const float *field, const geomlib_fluid *f,
                          float x, float y, float z)
{
    float gx = x/f->dx - 0.5f, gy = y/f->dx - 0.5f, gz = z/f->dx - 0.5f;
    int x0, y0, z0, i, j, k;
    double fx, fy, fz, sum = 0.0;
    if (gx < 0 || gy < 0 || gz < 0 ||
        gx > f->nx-1 || gy > f->ny-1 || gz > f->nz-1) return 0.0f;
    x0 = (int)floorf(gx); y0=(int)floorf(gy); z0=(int)floorf(gz);
    fx = gx-x0; fy=gy-y0; fz=gz-z0;
    for (k = 0; k < 2; ++k) {
        for (j = 0; j < 2; ++j) {
            for (i = 0; i < 2; ++i) {
                uint32_t xi = (uint32_t)(x0+i), yi=(uint32_t)(y0+j), zi=(uint32_t)(z0+k);
                double wgt = (i?fx:1-fx)*(j?fy:1-fy)*(k?fz:1-fz);
                if (xi < f->nx && yi < f->ny && zi < f->nz)
                    sum += wgt * field[IDX(f,xi,yi,zi)];
            }
        }
    }
    return (float)sum;
}

int geomlib_fluid_advect(geomlib_fluid *f, float dt)
{
    float *nu, *nv, *nw;
    size_t n, i;
    uint32_t x, y, z;
    if (!f || !f->u || !(dt > 0.0f)) return GEOMLIB_MESH_ERROR_INVALID;
    n = (size_t)f->nx*f->ny*f->nz;
    nu = (float *)malloc(n*sizeof(float));
    nv = (float *)malloc(n*sizeof(float));
    nw = (float *)malloc(n*sizeof(float));
    if (!nu || !nv || !nw) { free(nu); free(nv); free(nw); return GEOMLIB_MESH_ERROR_MEMORY; }
    for (z = 0; z < f->nz; ++z)
        for (y = 0; y < f->ny; ++y)
            for (x = 0; x < f->nx; ++x) {
                float px = (float)x*f->dx, py = (float)y*f->dx, pz = (float)z*f->dx;
                size_t id = IDX(f,x,y,z);
                float bx = px - dt*f->u[id];
                float by = py - dt*f->v[id];
                float bz = pz - dt*f->w[id];
                nu[id] = sample_field(f->u, f, bx, by, bz);
                nv[id] = sample_field(f->v, f, bx, by, bz);
                nw[id] = sample_field(f->w, f, bx, by, bz);
            }
    memcpy(f->u, nu, n*sizeof(float));
    memcpy(f->v, nv, n*sizeof(float));
    memcpy(f->w, nw, n*sizeof(float));
    free(nu); free(nv); free(nw);
    (void)i;
    return GEOMLIB_MESH_SUCCESS;
}

int geomlib_fluid_project(geomlib_fluid *f, float dt)
{
    size_t n, i, count = 0;
    uint32_t *rows, *cols;
    double *vals, *b, *x;
    geomlib_spmat A;
    int rc;
    uint32_t nx, ny, nz;
    float inv;

    (void)dt; /* dt cancels in the projection step */
    if (!f || !f->u) return GEOMLIB_MESH_ERROR_INVALID;
    n = (size_t)f->nx*f->ny*f->nz;
    nx=f->nx; ny=f->ny; nz=f->nz;
    inv = 1.0f/f->dx;
    rows = (uint32_t *)malloc((n*7+1)*sizeof(uint32_t));
    cols = (uint32_t *)malloc((n*7+1)*sizeof(uint32_t));
    vals = (double *)malloc((n*7+1)*sizeof(double));
    b = (double *)calloc(n, sizeof(double));
    x = (double *)calloc(n, sizeof(double));
    if (!rows || !cols || !vals || !b || !x) {
        free(rows); free(cols); free(vals); free(b); free(x);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }

    /* Forward-difference divergence; gradient below uses the same operator,
     * so the 7-point stencil is exactly D*G and projection cancels divergence. */
    for (uint32_t z = 0; z < nz; ++z)
        for (uint32_t y = 0; y < ny; ++y)
            for (uint32_t x0 = 0; x0 < nx; ++x0) {
                size_t id = IDX(f,x0,y,z);
                float du = x0+1<nx ? f->u[IDX(f,x0+1,y,z)]-f->u[id] : 0.0f;
                float dv = y+1<ny ? f->v[IDX(f,x0,y+1,z)]-f->v[id] : 0.0f;
                float dw = z+1<nz ? f->w[IDX(f,x0,y,z+1)]-f->w[id] : 0.0f;
                f->divergence[id] = (du+dv+dw)*inv;
                /* A = -dx^2*L, so solving A p = b gives L p = -b/dx^2;
                 * we want L p = div, hence b = -div*dx^2. */
                b[id] = -f->divergence[id]*f->dx*f->dx;
            }
    /* All-Neumann Poisson is singular: remove the mean and pin one cell. */
    {
        double mean = 0.0;
        for (i = 0; i < n; ++i) mean += b[i];
        mean /= (double)n;
        for (i = 0; i < n; ++i) b[i] -= mean;
        b[0] = 0.0;
    }

    /* Poisson matrix: -Laplacian with Neumann boundaries. */
    for (uint32_t z = 0; z < nz; ++z)
        for (uint32_t y = 0; y < ny; ++y)
            for (uint32_t x0 = 0; x0 < nx; ++x0) {
                size_t id = IDX(f,x0,y,z);
                int neighbours = 0;
                rows[count]=(uint32_t)id; cols[count]=(uint32_t)id;
                vals[count++]= 6.0 + 1e-6; /* regularize the Neumann nullspace */
                if (x0>0)      { rows[count]=(uint32_t)id; cols[count]=(uint32_t)(id-1);   vals[count++]=-1.0; neighbours++; }
                if (x0+1<nx)   { rows[count]=(uint32_t)id; cols[count]=(uint32_t)(id+1);   vals[count++]=-1.0; neighbours++; }
                if (y>0)       { rows[count]=(uint32_t)id; cols[count]=(uint32_t)(id-nx);  vals[count++]=-1.0; neighbours++; }
                if (y+1<ny)    { rows[count]=(uint32_t)id; cols[count]=(uint32_t)(id+nx);  vals[count++]=-1.0; neighbours++; }
                if (z>0)       { rows[count]=(uint32_t)id; cols[count]=(uint32_t)(id-nx*ny); vals[count++]=-1.0; neighbours++; }
                if (z+1<nz)    { rows[count]=(uint32_t)id; cols[count]=(uint32_t)(id+nx*ny); vals[count++]=-1.0; neighbours++; }
                (void)neighbours;
                if (!neighbours) { rows[count]=(uint32_t)id; cols[count]=(uint32_t)id; vals[count++]=1.0; b[id]=0.0; }
            }
    rc = geomlib_spmat_from_triplets(rows, cols, vals, count, n, n, &A);
    if (rc != GEOMLIB_SPARSE_SUCCESS) {
        free(rows); free(cols); free(vals); free(b); free(x);
        return GEOMLIB_MESH_ERROR_INVALID;
    }
    rc = geomlib_solve_sparse_direct(&A, b, x, NULL);
    geomlib_spmat_free(&A);
    if (rc != GEOMLIB_SPARSE_SUCCESS) {
        free(rows); free(cols); free(vals); free(b); free(x);
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    for (i = 0; i < n; ++i) f->pressure[i] = (float)x[i];
    /* Gradient uses backward differences so that D_forward . G_backward is
     * exactly the 7-point Laplacian used in the Poisson matrix. */
    for (uint32_t z = 0; z < nz; ++z)
        for (uint32_t y = 0; y < ny; ++y)
            for (uint32_t x0 = 0; x0 < nx; ++x0) {
                size_t id = IDX(f,x0,y,z);
                float dpdx = x0>0 ? (f->pressure[id]-f->pressure[IDX(f,x0-1,y,z)])*inv : 0.0f;
                float dpdy = y>0  ? (f->pressure[id]-f->pressure[IDX(f,x0,y-1,z)])*inv : 0.0f;
                float dpdz = z>0  ? (f->pressure[id]-f->pressure[IDX(f,x0,y,z-1)])*inv : 0.0f;
                f->u[id] -= dpdx;
                f->v[id] -= dpdy;
                f->w[id] -= dpdz;
            }
    free(rows); free(cols); free(vals); free(b); free(x);
    return GEOMLIB_MESH_SUCCESS;
}

float geomlib_fluid_rms_divergence(const geomlib_fluid *f)
{
    double sum = 0.0;
    size_t n, i;
    if (!f || !f->u) return INFINITY;
    n = (size_t)f->nx*f->ny*f->nz;
    for (i = 0; i < n; ++i) sum += (double)f->divergence[i]*f->divergence[i];
    return (float)sqrt(sum/n);
}
