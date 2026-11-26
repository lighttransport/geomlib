/*
 * geomlib - Feature-preserving point-set processing (see geomlib_pointset.h)
 *
 * Clean-room implementation of bilateral filtering for point clouds:
 *   x_i' = sum_j w_ij x_j / sum_j w_ij
 *   w_ij = exp(-||x_i-x_j||^2 / (2 sigma_d^2))
 *          * exp(-(1 - n_i . n_j) / sigma_n)
 * and of statistical outlier removal based on average k-NN distance.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_pointset.h"

#include <float.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_bilateral_params geomlib_bilateral_params_default(void)
{
    geomlib_bilateral_params p;
    p.spatial_sigma = 0.0f;
    p.normal_sigma = 0.5f;
    p.neighbors = 16;
    p.iterations = 2;
    return p;
}

static float spacing_estimate(const geomlib_pointcloud *cloud, int k)
{
    size_t i, j, count = 0;
    double sum = 0.0;
    if (cloud->point_count < 2) return 1e-3f;
    for (i = 0; i < cloud->point_count; i += (cloud->point_count > 512 ? 8 : 1)) {
        double best[2] = {DBL_MAX, DBL_MAX};
        for (j = 0; j < cloud->point_count; ++j) {
            geomlib_vec3 d = geomlib_vec3_sub(cloud->points[i], cloud->points[j]);
            double dd = geomlib_vec3_dot(d, d);
            if (i == j || dd < 1e-18) continue;
            if (dd < best[0]) { best[1] = best[0]; best[0] = dd; }
            else if (dd < best[1]) best[1] = dd;
        }
        if (best[1] < DBL_MAX) { sum += sqrt(best[0]) + sqrt(best[1]); count++; }
    }
    (void)k;
    return count ? (float)(sum / (2.0*count)) : 1e-3f;
}

int geomlib_pointcloud_bilateral_smooth(geomlib_pointcloud *cloud,
                                        const geomlib_bilateral_params *params)
{
    size_t n, i, j, iter;
    int k, iters;
    float sigma_d, sigma_n;
    geomlib_vec3 *source, *dest;
    int rc;

    if (!cloud || !cloud->points || cloud->point_count < 3)
        return GEOMLIB_MESH_ERROR_INVALID;
    n = cloud->point_count;
    k = params ? params->neighbors : 16;
    if (k < 2) k = 2;
    if (k > 64) k = 64;
    iters = params ? params->iterations : 2;
    if (iters < 0) iters = 0;
    sigma_d = params ? params->spatial_sigma : 0.0f;
    if (!(sigma_d > 0.0f)) sigma_d = spacing_estimate(cloud, k);
    if (!(sigma_d > 0.0f)) sigma_d = 1e-3f;
    sigma_n = params ? params->normal_sigma : 0.5f;
    if (!(sigma_n > 0.0f)) sigma_n = 0.5f;

    if (!cloud->normals) {
        rc = geomlib_pointcloud_estimate_normals(cloud, k);
        if (rc != GEOMLIB_MESH_SUCCESS) return rc;
    }

    source = (geomlib_vec3 *)malloc(n*sizeof(geomlib_vec3));
    dest = (geomlib_vec3 *)malloc(n*sizeof(geomlib_vec3));
    if (!source || !dest) { free(source); free(dest); return GEOMLIB_MESH_ERROR_MEMORY; }

    for (iter = 0; iter < (size_t)iters; ++iter) {
        memcpy(source, cloud->points, n*sizeof(geomlib_vec3));
        for (i = 0; i < n; ++i) {
            size_t *idx = (size_t *)malloc((size_t)k*sizeof(size_t));
            double *dist = (double *)malloc((size_t)k*sizeof(double));
            size_t count = 0;
            double wx=0, wy=0, wz=0, wsum=0;
            if (!idx || !dist) { free(idx); free(dist); continue; }
            for (j = 0; j < n; ++j) {
                geomlib_vec3 d;
                double dd;
                size_t pos;
                if (i == j) continue;
                d = geomlib_vec3_sub(source[i], source[j]);
                dd = geomlib_vec3_dot(d, d);
                if (count == (size_t)k && dd >= dist[count-1]) continue;
                pos = count < (size_t)k ? count : (size_t)(k-1);
                while (pos > 0 && dist[pos-1] > dd) {
                    dist[pos] = dist[pos-1];
                    idx[pos] = idx[pos-1];
                    --pos;
                }
                dist[pos] = dd;
                idx[pos] = j;
                if (count < (size_t)k) count++;
            }
            for (j = 0; j < count; ++j) {
                double ds = dist[j] / (2.0*sigma_d*sigma_d);
                double nd = 1.0 - geomlib_vec3_dot(cloud->normals[i],
                                                    cloud->normals[idx[j]]);
                double w = exp(-ds - nd/sigma_n);
                wx += w * source[idx[j]].x;
                wy += w * source[idx[j]].y;
                wz += w * source[idx[j]].z;
                wsum += w;
            }
            dest[i] = wsum > 0 ? geomlib_vec3_create((float)(wx/wsum),
                                                     (float)(wy/wsum), (float)(wz/wsum))
                               : source[i];
            free(idx); free(dist);
        }
        memcpy(cloud->points, dest, n*sizeof(geomlib_vec3));
    }
    free(source); free(dest);
    return GEOMLIB_MESH_SUCCESS;
}

typedef struct { size_t index; double dist; } rank_item;

static int cmp_rank(const void *a, const void *b)
{
    const rank_item *x = (const rank_item *)a, *y = (const rank_item *)b;
    return x->dist > y->dist ? -1 : (x->dist < y->dist ? 1 : 0);
}

size_t geomlib_pointcloud_remove_outliers(geomlib_pointcloud *cloud,
                                          float fraction, int k)
{
    size_t n, i, j, remove_count, kept = 0;
    rank_item *rank;
    bool *removed;

    if (!cloud || !cloud->points || cloud->point_count < 4 ||
        !(fraction > 0.0f) || fraction >= 1.0f)
        return 0;
    n = cloud->point_count;
    if (k < 2) k = 2;
    if (k > 64) k = 64;
    rank = (rank_item *)malloc(n*sizeof(rank_item));
    removed = (bool *)calloc(n, sizeof(bool));
    if (!rank || !removed) { free(rank); free(removed); return 0; }

    for (i = 0; i < n; ++i) {
        double best[4] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
        for (j = 0; j < n; ++j) {
            geomlib_vec3 d = geomlib_vec3_sub(cloud->points[i], cloud->points[j]);
            double dd = geomlib_vec3_dot(d, d);
            if (i == j) continue;
            if (dd < best[0]) {
                best[3]=best[2]; best[2]=best[1]; best[1]=best[0]; best[0]=dd;
            } else if (dd < best[1]) {
                best[3]=best[2]; best[2]=best[1]; best[1]=dd;
            } else if (dd < best[2]) {
                best[3]=best[2]; best[2]=dd;
            } else if (dd < best[3]) {
                best[3]=dd;
            }
        }
        rank[i].index = i;
        rank[i].dist = (best[0]+best[1]+best[2]+best[3]) / 4.0;
    }
    qsort(rank, n, sizeof(rank_item), cmp_rank);
    remove_count = (size_t)(n * fraction);
    if (remove_count >= n) remove_count = n-1;
    for (i = 0; i < remove_count; ++i) removed[rank[i].index] = true;

    for (i = 0; i < n; ++i) {
        if (removed[i]) continue;
        cloud->points[kept] = cloud->points[i];
        if (cloud->normals) cloud->normals[kept] = cloud->normals[i];
        if (cloud->confidences) cloud->confidences[kept] = cloud->confidences[i];
        kept++;
    }
    cloud->point_count = kept;
    free(rank); free(removed);
    return remove_count;
}
