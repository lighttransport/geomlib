/*
 * geomlib - DBSCAN (see geomlib_dbscan.h)
 *
 * Clean-room radius-neighbourhood expansion with a k-d tree.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_dbscan.h"
#include "geomlib/geomlib_kdtree.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

geomlib_dbscan_params geomlib_dbscan_params_default(void)
{
    geomlib_dbscan_params p;
    p.eps=0.0f;
    p.min_points=5;
    return p;
}

int geomlib_pointcloud_dbscan(const geomlib_pointcloud *cloud,
                              const geomlib_dbscan_params *params,
                              uint32_t *labels, size_t *cluster_count)
{
    geomlib_dbscan_params p=params ? *params : geomlib_dbscan_params_default();
    geomlib_kdtree *tree=NULL;
    uint32_t *nbrs=NULL;
    float *dst=NULL;
    uint32_t *queue=NULL;
    size_t n,i,k,clusters=0;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!cloud || !cloud->points || !cloud->point_count || !labels)
        return rc;
    n=cloud->point_count;
    if (p.min_points<2) p.min_points=2;
    for (i=0;i<n;++i) labels[i]=GEOMLIB_DBSCAN_NOISE;

    if (!(p.eps>0.0f)) {
        /* Estimate twice the mean nearest-neighbour distance. */
        double sum=0.0;
        tree=geomlib_kdtree_build(cloud->points,n);
        if (!tree) return GEOMLIB_MESH_ERROR_MEMORY;
        for (i=0;i<n;++i) {
            uint32_t id=0; float d=0.0f;
            if (geomlib_kdtree_knn(tree,cloud->points[i],2,&id,&d)>=2)
                sum+=d;
        }
        p.eps=(float)(2.0*sum/n);
        geomlib_kdtree_free(tree);
        tree=NULL;
    }
    if (!(p.eps>1e-12f)) return rc;

    tree=geomlib_kdtree_build(cloud->points,n);
    nbrs=(uint32_t*)malloc(n*sizeof(uint32_t));
    dst=(float*)malloc(n*sizeof(float));
    queue=(uint32_t*)malloc(n*sizeof(uint32_t));
    if (!tree||!nbrs||!dst||!queue) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto cleanup; }

    for (i=0;i<n;++i) {
        size_t found;
        if (labels[i]!=GEOMLIB_DBSCAN_NOISE) continue;
        found=geomlib_kdtree_radius(tree,cloud->points[i],p.eps,nbrs,dst,n);
        if (found<(size_t)p.min_points) continue;
        labels[i]=(uint32_t)clusters;
        {
            size_t head=0,tail=0;
            queue[tail++]=nbrs[0];
            for (k=1;k<found;++k) {
                if (nbrs[k]==i) continue;
                queue[tail++]=nbrs[k];
            }
            while (head<tail) {
                uint32_t q=queue[head++];
                size_t cnt,f;
                if (labels[q]==GEOMLIB_DBSCAN_NOISE) labels[q]=(uint32_t)clusters;
                if (labels[q]!=(uint32_t)clusters) continue;
                cnt=geomlib_kdtree_radius(tree,cloud->points[q],p.eps,nbrs,dst,n);
                if (cnt<(size_t)p.min_points) continue;
                for (f=0;f<cnt;++f) {
                    uint32_t r=nbrs[f];
                    if (labels[r]==GEOMLIB_DBSCAN_NOISE) {
                        labels[r]=(uint32_t)clusters;
                        if (tail<n) queue[tail++]=r;
                    }
                }
            }
        }
        clusters++;
    }
    if (cluster_count) *cluster_count=clusters;
    rc=GEOMLIB_MESH_SUCCESS;

cleanup:
    geomlib_kdtree_free(tree);
    free(nbrs); free(dst); free(queue);
    return rc;
}
