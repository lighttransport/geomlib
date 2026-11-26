/*
 * geomlib - Point-cloud k-means (see geomlib_kmeans.h)
 *
 * Clean-room deterministic farthest-point seeding and Lloyd iterations.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_kmeans.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double dist2(geomlib_vec3 a, geomlib_vec3 b)
{
    double dx=(double)a.x-b.x, dy=(double)a.y-b.y, dz=(double)a.z-b.z;
    return dx*dx+dy*dy+dz*dz;
}

int geomlib_pointcloud_kmeans(const geomlib_pointcloud *cloud, size_t k,
                              int max_iterations,
                              uint32_t *labels, geomlib_vec3 *centroids)
{
    geomlib_vec3 *c=NULL;
    size_t i,j,n;
    int iter;

    if (!cloud || !cloud->points || !cloud->point_count || !labels || !k)
        return GEOMLIB_MESH_ERROR_INVALID;
    n=cloud->point_count;
    if (k>n) k=n;
    if (max_iterations<1) max_iterations=1;
    for (i=0;i<n;++i) labels[i]=0;
    c=(geomlib_vec3*)malloc(k*sizeof(geomlib_vec3));
    if (!c) return GEOMLIB_MESH_ERROR_MEMORY;

    /* Seed at the point farthest from the current centroids. */
    c[0]=cloud->points[0];
    for (j=1;j<k;++j) {
        size_t far=0;
        double bd=-1.0;
        for (i=0;i<n;++i) {
            double md=DBL_MAX;
            for (size_t q=0;q<j;++q) {
                double d=dist2(cloud->points[i],c[q]);
                if (d<md) md=d;
            }
            if (md>bd) { bd=md; far=i; }
        }
        c[j]=cloud->points[far];
    }

    for (iter=0;iter<max_iterations;++iter) {
        double *sx=calloc(k,sizeof(double));
        double *sy=calloc(k,sizeof(double));
        double *sz=calloc(k,sizeof(double));
        size_t *count=calloc(k,sizeof(size_t));
        int changed=0;
        if (!sx||!sy||!sz||!count) {
            free(sx);free(sy);free(sz);free(count);free(c);
            return GEOMLIB_MESH_ERROR_MEMORY;
        }
        for (i=0;i<n;++i) {
            size_t best=0;
            double bd=DBL_MAX;
            for (j=0;j<k;++j) {
                double d=dist2(cloud->points[i],c[j]);
                if (d<bd) { bd=d; best=j; }
            }
            if (labels[i]!=best) { labels[i]=(uint32_t)best; changed=1; }
            sx[best]+=cloud->points[i].x;
            sy[best]+=cloud->points[i].y;
            sz[best]+=cloud->points[i].z;
            count[best]++;
        }
        for (j=0;j<k;++j)
            if (count[j])
                c[j]=geomlib_vec3_create((float)(sx[j]/count[j]),
                                         (float)(sy[j]/count[j]),
                                         (float)(sz[j]/count[j]));
        free(sx);free(sy);free(sz);free(count);
        if (!changed) break;
    }
    if (centroids) memcpy(centroids,c,k*sizeof(geomlib_vec3));
    free(c);
    return (int)k;
}
