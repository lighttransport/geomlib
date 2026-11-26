/*
 * geomlib - Poisson surface reconstruction (see geomlib_poisson.h)
 *
 * Voxel-based clean-room implementation:
 *   1. Splat oriented normals into a vector field on a regular grid.
 *   2. Compute the discrete divergence of that field.
 *   3. Solve the discrete Poisson equation L * chi = h^2 * div(chi) with a
 *      Jacobi-preconditioned CG solver (one corner pinned to remove the
 *      constant nullspace).
 *   4. Choose the iso level as the median indicator value at the input points
 *      and extract the surface with the existing marching-cubes path.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_poisson.h"
#include "geomlib/geomlib_sdf.h"
#include "geomlib/geomlib_sparse.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_poisson_params geomlib_poisson_params_default(void)
{
    geomlib_poisson_params p;
    p.resolution = 48;
    p.padding = 0.1f;
    p.max_iterations = 1000;
    p.tolerance = 1e-6f;
    return p;
}

static size_t idx3(size_t x, size_t y, size_t z, size_t nx, size_t ny)
{
    return x + nx * (y + ny * z);
}

/* Trilinear sample of a scalar grid at a world-space point. */
static double sample_grid(const double *values, const geomlib_sdf *grid, geomlib_vec3 p)
{
    double gx = ((double)p.x - grid->origin.x) / grid->voxel_size - 0.5;
    double gy = ((double)p.y - grid->origin.y) / grid->voxel_size - 0.5;
    double gz = ((double)p.z - grid->origin.z) / grid->voxel_size - 0.5;
    size_t nx = grid->resolution[0], ny = grid->resolution[1], nz = grid->resolution[2];
    int x0, y0, z0;
    double fx, fy, fz, sum = 0.0;
    int i, j, k;

    if (gx < 0 || gy < 0 || gz < 0 ||
        gx > (double)nx-1 || gy > (double)ny-1 || gz > (double)nz-1)
        return 0.0;
    x0 = (int)floor(gx); y0 = (int)floor(gy); z0 = (int)floor(gz);
    fx = gx - x0; fy = gy - y0; fz = gz - z0;
    for (k = 0; k < 2; ++k) {
        for (j = 0; j < 2; ++j) {
            for (i = 0; i < 2; ++i) {
                size_t xi = (size_t)(x0+i), yi = (size_t)(y0+j), zi = (size_t)(z0+k);
                double w = (i ? fx : 1.0-fx) * (j ? fy : 1.0-fy) * (k ? fz : 1.0-fz);
                if (xi < nx && yi < ny && zi < nz)
                    sum += w * values[idx3(xi, yi, zi, nx, ny)];
            }
        }
    }
    return sum;
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

int geomlib_reconstruct_poisson(const geomlib_pointcloud *cloud,
                                const geomlib_poisson_params *params,
                                geomlib_trimesh *out)
{
    geomlib_aabb bounds = {0};
    geomlib_vec3 size, center;
    double max_extent;
    float h;
    uint32_t res[3];
    size_t nx, ny, nz, total, i;
    float *fx = NULL, *fy = NULL, *fz = NULL;
    double *div = NULL, *chi = NULL, *at_points = NULL;
    uint32_t *rows = NULL, *cols = NULL;
    double *vals = NULL;
    geomlib_spmat l;
    geomlib_sparse_solve_params sp;
    geomlib_sdf grid;
    double iso;
    int rc = GEOMLIB_MESH_ERROR_INVALID;
    uint32_t resmax;

    if (!cloud || !out || !cloud->points || cloud->point_count == 0)
        return GEOMLIB_MESH_ERROR_INVALID;
    if (!cloud->normals) return GEOMLIB_MESH_ERROR_INVALID;

    /* Compute point bounds directly. */
    {
        geomlib_vec3 lo = cloud->points[0], hi = cloud->points[0];
        for (i = 1; i < cloud->point_count; ++i) {
            geomlib_vec3 p = cloud->points[i];
            if (p.x < lo.x) lo.x = p.x;
            if (p.y < lo.y) lo.y = p.y;
            if (p.z < lo.z) lo.z = p.z;
            if (p.x > hi.x) hi.x = p.x;
            if (p.y > hi.y) hi.y = p.y;
            if (p.z > hi.z) hi.z = p.z;
        }
        bounds.min = lo;
        bounds.max = hi;
    }
    size = geomlib_aabb_size(bounds);
    center = geomlib_aabb_center(bounds);
    max_extent = fmax(size.x, fmax(size.y, size.z));
    if (!(max_extent > 0.0)) return GEOMLIB_MESH_ERROR_INVALID;

    resmax = params ? params->resolution : 48;
    if (resmax < 16) resmax = 16;
    if (resmax > 128) resmax = 128;
    h = (float)(max_extent / resmax);
    {
        float pad = params ? params->padding : 0.1f;
        if (pad < 0.0f) pad = 0.0f;
        if (pad > 1.0f) pad = 1.0f;
        res[0] = (uint32_t)(size.x / h) + 5 + (uint32_t)(2.0f * pad * size.x / h);
        res[1] = (uint32_t)(size.y / h) + 5 + (uint32_t)(2.0f * pad * size.y / h);
        res[2] = (uint32_t)(size.z / h) + 5 + (uint32_t)(2.0f * pad * size.z / h);
    }
    if (res[0] > 128) res[0] = 128;
    if (res[1] > 128) res[1] = 128;
    if (res[2] > 128) res[2] = 128;
    nx = res[0]; ny = res[1]; nz = res[2];
    total = nx * ny * nz;

    fx = (float *)calloc(total, sizeof(float));
    fy = (float *)calloc(total, sizeof(float));
    fz = (float *)calloc(total, sizeof(float));
    div = (double *)calloc(total, sizeof(double));
    chi = (double *)calloc(total, sizeof(double));
    at_points = (double *)malloc(cloud->point_count * sizeof(double));
    rows = (uint32_t *)malloc((total * 7 + 1) * sizeof(uint32_t));
    cols = (uint32_t *)malloc((total * 7 + 1) * sizeof(uint32_t));
    vals = (double *)malloc((total * 7 + 1) * sizeof(double));
    if (!fx || !fy || !fz || !div || !chi || !at_points || !rows || !cols || !vals) {
        rc = GEOMLIB_MESH_ERROR_MEMORY;
        goto cleanup;
    }

    memset(&grid, 0, sizeof(grid));
    grid.voxel_size = h;
    grid.resolution[0] = (uint32_t)nx;
    grid.resolution[1] = (uint32_t)ny;
    grid.resolution[2] = (uint32_t)nz;
    grid.origin = geomlib_vec3_sub(center, geomlib_vec3_mul(size, 0.5f + (params ? params->padding : 0.1f)));
    grid.extent = (geomlib_vec3){ h*nx, h*ny, h*nz };
    grid.bounds.min = grid.origin;
    grid.bounds.max = geomlib_vec3_add(grid.origin, grid.extent);

    /* Splat oriented normals with trilinear weights (radius 2 voxels). */
    for (i = 0; i < cloud->point_count; ++i) {
        geomlib_vec3 p = cloud->points[i];
        geomlib_vec3 n = cloud->normals[i];
        double gx = ((double)p.x - grid.origin.x)/h - 0.5;
        double gy = ((double)p.y - grid.origin.y)/h - 0.5;
        double gz = ((double)p.z - grid.origin.z)/h - 0.5;
        int x0 = (int)floor(gx), y0 = (int)floor(gy), z0 = (int)floor(gz);
        int a, b, c;
        for (c = -1; c <= 2; ++c) {
            for (b = -1; b <= 2; ++b) {
                for (a = -1; a <= 2; ++a) {
                    int xi = x0+a, yi = y0+b, zi = z0+c;
                    if (xi < 0 || yi < 0 || zi < 0 ||
                        (size_t)xi >= nx || (size_t)yi >= ny || (size_t)zi >= nz) continue;
                    double dx = gx - xi, dy = gy - yi, dz = gz - zi;
                    double w = exp(-(dx*dx+dy*dy+dz*dz) * 0.5);
                    size_t k = idx3((size_t)xi, (size_t)yi, (size_t)zi, nx, ny);
                    fx[k] += (float)(w * n.x);
                    fy[k] += (float)(w * n.y);
                    fz[k] += (float)(w * n.z);
                }
            }
        }
    }

    /* Divergence with clamped central differences. */
    {
        size_t x, y, z;
        double inv2h = 1.0 / (2.0 * h);
        for (z = 0; z < nz; ++z) {
            for (y = 0; y < ny; ++y) {
                for (x = 0; x < nx; ++x) {
                    size_t k = idx3(x,y,z,nx,ny);
                    double dx = (fx[idx3(x+1<nx?x+1:x, y, z, nx, ny)] -
                                 fx[idx3(x>0?x-1:0, y, z, nx, ny)]) *
                                ((x==0 || x+1==nx) ? 1.0 : 0.5);
                    double dy = (fy[idx3(x, y+1<ny?y+1:y, z, nx, ny)] -
                                 fy[idx3(x, y>0?y-1:0, z, nx, ny)]) *
                                ((y==0 || y+1==ny) ? 1.0 : 0.5);
                    double dz = (fz[idx3(x, y, z+1<nz?z+1:z, nx, ny)] -
                                 fz[idx3(x, y, z>0?z-1:0, nx, ny)]) *
                                ((z==0 || z+1==nz) ? 1.0 : 0.5);
                    div[k] = (dx + dy + dz) * inv2h;
                }
            }
        }
    }

    /* Assemble the 7-point Laplacian with Neumann boundary conditions. */
    {
        size_t count = 0, x, y, z;
        for (z = 0; z < nz; ++z) {
            for (y = 0; y < ny; ++y) {
                for (x = 0; x < nx; ++x) {
                    size_t k = idx3(x,y,z,nx,ny);
                    int neighbors = 0;
                    rows[count] = (uint32_t)k;
                    cols[count] = (uint32_t)k;
                    vals[count++] = -6.0;
                    if (x > 0) { rows[count]=(uint32_t)k; cols[count]=(uint32_t)(k-1); vals[count++]=1.0; neighbors++; }
                    if (x+1<nx) { rows[count]=(uint32_t)k; cols[count]=(uint32_t)(k+1); vals[count++]=1.0; neighbors++; }
                    if (y > 0) { rows[count]=(uint32_t)k; cols[count]=(uint32_t)(k-nx); vals[count++]=1.0; neighbors++; }
                    if (y+1<ny) { rows[count]=(uint32_t)k; cols[count]=(uint32_t)(k+nx); vals[count++]=1.0; neighbors++; }
                    if (z > 0) { rows[count]=(uint32_t)k; cols[count]=(uint32_t)(k-nx*ny); vals[count++]=1.0; neighbors++; }
                    if (z+1<nz) { rows[count]=(uint32_t)k; cols[count]=(uint32_t)(k+nx*ny); vals[count++]=1.0; neighbors++; }
                    (void)neighbors;
                    /* Fix the corner value to zero to remove the nullspace. */
                    if (x == 0 && y == 0 && z == 0) {
                        rows[count]=(uint32_t)k; cols[count]=(uint32_t)k; vals[count++]=1.0;
                    }
                }
            }
        }
        rc = geomlib_spmat_from_triplets(rows, cols, vals, count, total, total, &l);
        if (rc != GEOMLIB_SPARSE_SUCCESS) { rc = GEOMLIB_MESH_ERROR_INVALID; goto cleanup; }
    }

    sp = geomlib_sparse_solve_params_default();
    {
        double *b = (double *)malloc(total * sizeof(double));
        if (!b) { geomlib_spmat_free(&l); rc = GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }
        for (i = 0; i < total; ++i) b[i] = div[i] * h * h;
        b[0] = 0.0;
        sp.max_iterations = params ? params->max_iterations : 1000;
        sp.tolerance = params ? params->tolerance : 1e-6f;
        rc = geomlib_solve_cg_sparse(&l, b, chi, &sp);
        free(b);
        geomlib_spmat_free(&l);
        if (rc != GEOMLIB_SPARSE_SUCCESS) { rc = GEOMLIB_MESH_ERROR_INVALID; goto cleanup; }
    }

    /* Iso level: median indicator value at the input points. */
    for (i = 0; i < cloud->point_count; ++i)
        at_points[i] = sample_grid(chi, &grid, cloud->points[i]);
    qsort(at_points, cloud->point_count, sizeof(double), cmp_double);
    iso = at_points[cloud->point_count / 2];

    grid.values = (float *)malloc(total * sizeof(float));
    if (!grid.values) { rc = GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }
    grid.total_voxels = total;
    for (i = 0; i < total; ++i) grid.values[i] = (float)(chi[i] - iso);
    rc = geomlib_sdf_extract_surface(&grid, out);
    free(grid.values);

cleanup:
    free(fx); free(fy); free(fz); free(div); free(chi);
    free(at_points); free(rows); free(cols); free(vals);
    return rc;
}
