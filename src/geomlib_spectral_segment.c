/*
 * geomlib - Spectral segmentation (see geomlib_spectral_segment.h)
 *
 * Clean-room harmonics embedding + farthest-seeded k-means.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_spectral_segment.h"
#include "geomlib/geomlib_harmonics.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_spectral_segment(const geomlib_trimesh *mesh, size_t k,
                             size_t eigen_count, uint32_t *labels)
{
    geomlib_harmonics_result h;
    double *centroids;
    size_t *counts;
    size_t n, dim, i, j, d, iter;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!mesh || !mesh->vertices || !mesh->triangle_count || !labels ||
        k<2 || eigen_count<2) return rc;
    n=mesh->vertex_count;
    if (k>n) k=n;
    if (eigen_count<k) eigen_count=k;
    dim=eigen_count-1;  /* skip the constant eigenvector */

    h=geomlib_harmonics_result_create(eigen_count,n);
    if (!h.eigenvalues || !h.eigenvectors) { geomlib_harmonics_result_free(&h);
        return GEOMLIB_MESH_ERROR_MEMORY; }
    rc=geomlib_compute_manifold_harmonics(mesh,eigen_count,&h);
    if (rc!=GEOMLIB_MESH_SUCCESS) { geomlib_harmonics_result_free(&h); return rc; }

    centroids=(double*)calloc(k*dim,sizeof(double));
    counts=(size_t*)calloc(k,sizeof(size_t));
    if (!centroids||!counts) { rc=GEOMLIB_MESH_ERROR_MEMORY; goto done; }

    /* Farthest-point seeding in embedding space. */
    for (i=0;i<k;++i) {
        size_t best=0;
        double bestd=-1.0;
        for (j=0;j<n;++j) {
            double dd=0.0;
            for (d=0;d<dim;++d) {
                double v=h.eigenvectors[(d+1)*n+j];
                dd+=v*v;
            }
            if (dd>bestd) { bestd=dd; best=j; }
        }
        for (d=0;d<dim;++d) centroids[i*dim+d]=h.eigenvectors[(d+1)*n+best];
        /* Avoid duplicate seeds by nudging along an axis. */
        for (j=0;j<i;++j) {
            double diff=0.0;
            for (d=0;d<dim;++d)
                diff+=fabs(centroids[i*dim+d]-centroids[j*dim+d]);
            if (diff<1e-12) centroids[i*dim+0]+=1e-6*(double)(i+1);
        }
    }

    for (iter=0;iter<32;++iter) {
        int changed=0;
        for (i=0;i<n;++i) {
            double bd=DBL_MAX;
            size_t bc=0;
            for (j=0;j<k;++j) {
                double dd=0.0;
                for (d=0;d<dim;++d) {
                    double v=h.eigenvectors[(d+1)*n+i]-centroids[j*dim+d];
                    dd+=v*v;
                }
                if (dd<bd) { bd=dd; bc=j; }
            }
            if (labels[i]!=bc) { labels[i]=(uint32_t)bc; changed=1; }
        }
        memset(counts,0,k*sizeof(size_t));
        memset(centroids,0,k*dim*sizeof(double));
        for (i=0;i<n;++i) {
            size_t c=labels[i];
            counts[c]++;
            for (d=0;d<dim;++d) centroids[c*dim+d]+=h.eigenvectors[(d+1)*n+i];
        }
        for (j=0;j<k;++j)
            if (counts[j])
                for (d=0;d<dim;++d) centroids[j*dim+d]/=(double)counts[j];
        if (!changed) break;
    }
    rc=GEOMLIB_MESH_SUCCESS;

done:
    free(centroids); free(counts);
    geomlib_harmonics_result_free(&h);
    return rc;
}
