/*
 * geomlib - RBF variational reconstruction (see geomlib_rbf_reconstruct.h)
 *
 * Clean-room implementation using compactly supported Wendland C2 kernels:
 *   phi(r) = (1-r)^4 (4r+1), r < 1; 0 otherwise.
 * Weights satisfy sum_j w_j phi(||p_i - p_j|| / support) = 1, so the implicit
 * field is positive inside the point cloud and crosses zero near the surface.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_rbf_reconstruct.h"
#include "geomlib/geomlib_sdf.h"
#include "geomlib/geomlib_solver.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_rbf_params geomlib_rbf_params_default(void)
{
    geomlib_rbf_params p;
    p.resolution = 40;
    p.support = 0.0f;
    p.max_points = 400;
    return p;
}

static float phi(float r)
{
    if (!(r < 1.0f)) return 0.0f;
    {
        float t = 1.0f-r;
        return t*t*t*t*(4.0f*r+1.0f);
    }
}

int geomlib_reconstruct_rbf(const geomlib_pointcloud *cloud,
                            const geomlib_rbf_params *params,
                            geomlib_trimesh *out)
{
    geomlib_vec3 *p;
    double *A, *w;
    float *values;
    geomlib_sdf grid;
    geomlib_aabb box;
    geomlib_vec3 size, center;
    size_t n, i, j, x, y, z;
    float support, spacing, diag;
    uint32_t res;
    int rc;

    if (!cloud || !out || !cloud->points || cloud->point_count < 4)
        return GEOMLIB_MESH_ERROR_INVALID;
    n = cloud->point_count;
    if (params && params->max_points > 0 && n > (size_t)params->max_points)
        n = (size_t)params->max_points;
    if (n > 500) n = 500;

    box = geomlib_trimesh_compute_aabb(NULL); /* replaced below */
    {
        geomlib_vec3 lo=cloud->points[0], hi=cloud->points[0];
        for (i = 1; i < n; ++i) {
            if (cloud->points[i].x<lo.x) lo.x=cloud->points[i].x;
            if (cloud->points[i].y<lo.y) lo.y=cloud->points[i].y;
            if (cloud->points[i].z<lo.z) lo.z=cloud->points[i].z;
            if (cloud->points[i].x>hi.x) hi.x=cloud->points[i].x;
            if (cloud->points[i].y>hi.y) hi.y=cloud->points[i].y;
            if (cloud->points[i].z>hi.z) hi.z=cloud->points[i].z;
        }
        box = geomlib_aabb_create(lo, hi);
    }
    size = geomlib_aabb_size(box);
    diag = geomlib_vec3_length(size);
    if (!(diag > 1e-9f)) return GEOMLIB_MESH_ERROR_INVALID;
    center = geomlib_aabb_center(box);
    spacing = diag / (float)cbrt((double)n);

    support = params ? params->support : 0.0f;
    if (!(support > 0.0f)) support = 3.0f*spacing;
    res = params ? params->resolution : 40;
    if (res < 16) res = 16;
    if (res > 96) res = 96;

    p = (geomlib_vec3 *)malloc(n*sizeof(geomlib_vec3));
    A = (double *)malloc(n*n*sizeof(double));
    w = (double *)malloc(n*sizeof(double));
    values = (float *)malloc((size_t)res*res*res*sizeof(float));
    if (!p || !A || !w || !values) {
        free(p); free(A); free(w); free(values);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    for (i = 0; i < n; ++i) p[i] = cloud->points[i];
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            geomlib_vec3 d = geomlib_vec3_sub(p[i], p[j]);
            A[i*n+j] = phi(geomlib_vec3_length(d)/support);
        }
        w[i] = 1.0;
    }
    rc = geomlib_solve_gauss_d(n, A, w);
    if (rc != GEOMLIB_SOLVER_SUCCESS) {
        free(p); free(A); free(w); free(values);
        return GEOMLIB_MESH_ERROR_INVALID;
    }

    /* Build a dense grid spanning the points with a small margin. */
    memset(&grid, 0, sizeof(grid));
    grid.voxel_size = 1.2f*diag/(float)res;
    grid.origin = geomlib_vec3_sub(center,
        geomlib_vec3_mul(size, 0.6f));
    grid.resolution[0] = grid.resolution[1] = grid.resolution[2] = res;
    grid.total_voxels = (size_t)res*res*res;
    grid.extent = geomlib_vec3_create(grid.voxel_size*res,
                                      grid.voxel_size*res, grid.voxel_size*res);
    grid.bounds.min = grid.origin;
    grid.bounds.max = geomlib_vec3_add(grid.origin, grid.extent);

    for (z = 0; z < res; ++z) {
        for (y = 0; y < res; ++y) {
            for (x = 0; x < res; ++x) {
                geomlib_vec3 q = {
                    grid.origin.x + ((float)x+0.5f)*grid.voxel_size,
                    grid.origin.y + ((float)y+0.5f)*grid.voxel_size,
                    grid.origin.z + ((float)z+0.5f)*grid.voxel_size
                };
                double f = 0.0;
                for (j = 0; j < n; ++j) {
                    geomlib_vec3 d = geomlib_vec3_sub(q, p[j]);
                    f += w[j]*phi(geomlib_vec3_length(d)/support);
                }
                values[x + res*(y + res*z)] = (float)(0.95 - f);
            }
        }
    }
    grid.values = values;
    rc = geomlib_sdf_extract_surface(&grid, out);
    free(grid.values);
    free(p); free(A); free(w);
    return rc;
}
