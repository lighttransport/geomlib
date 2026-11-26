/*
 * geomlib - MLS surface projection (see geomlib_mls.h)
 *
 * Clean-room implementation: weighted local PCA (cyclic Jacobi eigen solve)
 * and iterative projection onto the fitted plane.
 *
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_mls.h"
#include "geomlib/geomlib_kdtree.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_mls_params geomlib_mls_params_default(void)
{
    geomlib_mls_params p;
    p.neighbors = 16;
    p.iterations = 3;
    p.sigma = 0.0f;
    return p;
}

/* Smallest eigenvector of a symmetric 3x3 matrix by cyclic Jacobi rotations. */
static void smallest_eigenvector(double a[3][3], double out[3])
{
    double v[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    int sweep, p, q;
    for (sweep = 0; sweep < 32; ++sweep) {
        double off = fabs(a[0][1])+fabs(a[0][2])+fabs(a[1][2]);
        if (off < 1e-18) break;
        for (p = 0; p < 2; ++p)
            for (q = p+1; q < 3; ++q) {
                double apq = a[p][q];
                if (fabs(apq) < 1e-20) continue;
                double tau = (a[q][q]-a[p][p])/(2.0*apq);
                double t = (tau >= 0.0 ? 1.0 : -1.0) /
                           (fabs(tau)+sqrt(1.0+tau*tau));
                double c = 1.0/sqrt(1.0+t*t), s = t*c;
                double app=a[p][p], aqq=a[q][q], apq2=a[p][q];
                int k;
                a[p][p] = c*c*app - 2.0*s*c*apq2 + s*s*aqq;
                a[q][q] = s*s*app + 2.0*s*c*apq2 + c*c*aqq;
                a[p][q] = a[q][p] = 0.0;
                for (k = 0; k < 3; ++k) {
                    if (k!=p && k!=q) {
                        double akp=a[k][p], akq=a[k][q];
                        a[k][p]=a[p][k]=c*akp-s*akq;
                        a[k][q]=a[q][k]=s*akp+c*akq;
                    }
                    {
                        double vkp=v[k][p], vkq=v[k][q];
                        v[k][p]=c*vkp-s*vkq;
                        v[k][q]=s*vkp+c*vkq;
                    }
                }
            }
    }
    {
        int min = 0;
        if (a[1][1] < a[min][min]) min = 1;
        if (a[2][2] < a[min][min]) min = 2;
        out[0]=v[0][min]; out[1]=v[1][min]; out[2]=v[2][min];
    }
    {
        double len = sqrt(out[0]*out[0]+out[1]*out[1]+out[2]*out[2]);
        if (len > 0) { out[0]/=len; out[1]/=len; out[2]/=len; }
        else { out[0]=0; out[1]=0; out[2]=1; }
    }
}

int geomlib_mls_project(const geomlib_pointcloud *cloud, geomlib_vec3 query,
                        const geomlib_mls_params *params,
                        geomlib_vec3 *projected, geomlib_vec3 *normal)
{
    geomlib_mls_params p = params ? *params : geomlib_mls_params_default();
    geomlib_kdtree *tree;
    uint32_t *idx;
    float *dst;
    geomlib_vec3 q = query, fitted_normal = geomlib_vec3_create(0,0,1);
    size_t k, iter;
    int rc = GEOMLIB_MESH_ERROR_INVALID;

    if (!cloud || !cloud->points || !cloud->point_count ||
        !projected || p.neighbors < 3) return rc;
    k = (size_t)p.neighbors;
    if (k > cloud->point_count) k = cloud->point_count;
    if (p.iterations < 1) p.iterations = 1;
    if (p.iterations > 16) p.iterations = 16;

    tree = geomlib_kdtree_build(cloud->points, cloud->point_count);
    idx = (uint32_t *)malloc(k*sizeof(uint32_t));
    dst = (float *)malloc(k*sizeof(float));
    if (!tree || !idx || !dst) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }

    for (iter = 0; iter < (size_t)p.iterations; ++iter) {
        size_t found, i;
        double centroid[3]={0,0,0}, wsum=0.0, sigma;
        double cov[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        double ev[3];
        found = geomlib_kdtree_knn(tree, q, k, idx, dst);
        if (found < 3) { rc=GEOMLIB_MESH_ERROR_INVALID; goto cleanup; }
        sigma = p.sigma > 0.0f ? p.sigma : (dst[found-1] > 1e-6f ? dst[found-1] : 1e-3);
        for (i = 0; i < found; ++i) {
            geomlib_vec3 v = cloud->points[idx[i]];
            double d = (double)dst[i]/sigma;
            double w = exp(-d*d);
            centroid[0]+=w*v.x; centroid[1]+=w*v.y; centroid[2]+=w*v.z;
            wsum+=w;
        }
        if (!(wsum > 1e-18)) { rc=GEOMLIB_MESH_ERROR_INVALID; goto cleanup; }
        centroid[0]/=wsum; centroid[1]/=wsum; centroid[2]/=wsum;
        for (i = 0; i < found; ++i) {
            geomlib_vec3 v = cloud->points[idx[i]];
            double dx=(double)v.x-centroid[0], dy=(double)v.y-centroid[1],
                   dz=(double)v.z-centroid[2];
            double d=(double)dst[i]/sigma, w=exp(-d*d);
            cov[0][0]+=w*dx*dx; cov[0][1]+=w*dx*dy; cov[0][2]+=w*dx*dz;
            cov[1][1]+=w*dy*dy; cov[1][2]+=w*dy*dz; cov[2][2]+=w*dz*dz;
        }
        cov[1][0]=cov[0][1]; cov[2][0]=cov[0][2]; cov[2][1]=cov[1][2];
        smallest_eigenvector(cov, ev);
        fitted_normal = geomlib_vec3_create((float)ev[0],(float)ev[1],(float)ev[2]);
        {
            double dot = ((double)q.x-centroid[0])*ev[0]
                       + ((double)q.y-centroid[1])*ev[1]
                       + ((double)q.z-centroid[2])*ev[2];
            q = geomlib_vec3_sub(q, geomlib_vec3_mul(fitted_normal,(float)dot));
        }
    }
    *projected = q;
    if (normal) *normal = fitted_normal;
    rc = GEOMLIB_MESH_SUCCESS;

cleanup:
    geomlib_kdtree_free(tree);
    free(idx); free(dst);
    return rc;
}
