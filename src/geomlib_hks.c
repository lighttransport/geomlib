/*
 * geomlib - Heat kernel signature (see geomlib_hks.h)
 *
 * Clean-room implementation from the Laplace-Beltrami spectrum.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_hks.h"
#include "geomlib/geomlib_harmonics.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_compute_hks(const geomlib_trimesh *mesh, size_t eigen_count,
                        const double *times, size_t ntimes, double *out)
{
    geomlib_harmonics_result h;
    size_t i, k, ti;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!mesh || !mesh->vertices || !mesh->triangle_count ||
        !times || !out || !ntimes || !eigen_count) return rc;
    for (ti=0;ti<ntimes;++ti)
        if (!(times[ti]>0.0) || !isfinite(times[ti])) return rc;

    h=geomlib_harmonics_result_create(eigen_count,mesh->vertex_count);
    if (!h.eigenvalues || !h.eigenvectors) {
        geomlib_harmonics_result_free(&h);
        return GEOMLIB_MESH_ERROR_MEMORY;
    }
    rc=geomlib_compute_manifold_harmonics(mesh,eigen_count,&h);
    if (rc!=GEOMLIB_MESH_SUCCESS) {
        geomlib_harmonics_result_free(&h);
        return rc;
    }

    for (ti=0;ti<ntimes;++ti) {
        double mx=0.0;
        for (i=0;i<mesh->vertex_count;++i) {
            double sum=0.0;
            for (k=0;k<h.count;++k) {
                double phi=h.eigenvectors[k*h.vertex_count+i];
                double e=exp(-h.eigenvalues[k]*times[ti]);
                sum+=e*phi*phi;
            }
            out[i*ntimes+ti]=sum;
            if (sum>mx) mx=sum;
        }
        /* Normalize each time slice by its maximum, as is standard for HKS. */
        if (mx>0.0)
            for (i=0;i<mesh->vertex_count;++i) out[i*ntimes+ti]/=mx;
    }
    geomlib_harmonics_result_free(&h);
    return GEOMLIB_MESH_SUCCESS;
}
