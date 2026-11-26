/*
 * geomlib - Discrete Fréchet distance (see geomlib_frechet.h)
 *
 * Clean-room dynamic programming. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_frechet.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double geomlib_frechet_distance(const geomlib_vec3 *a, size_t na,
                                const geomlib_vec3 *b, size_t nb)
{
    double *prev=NULL,*cur=NULL,result;
    size_t i,j;

    if (!a||!b||!na||!nb) return INFINITY;
    prev=(double*)malloc(nb*sizeof(double));
    cur=(double*)malloc(nb*sizeof(double));
    if (!prev||!cur) { free(prev); free(cur); return INFINITY; }
    for (j=0;j<nb;++j)
        prev[j]=geomlib_vec3_length(geomlib_vec3_sub(a[0],b[j]));
    for (i=1;i<na;++i) {
        cur[0]=fmax(prev[0],geomlib_vec3_length(geomlib_vec3_sub(a[i],b[0])));
        for (j=1;j<nb;++j) {
            double d=geomlib_vec3_length(geomlib_vec3_sub(a[i],b[j]));
            double m=fmin(fmax(cur[j-1],prev[j]),fmax(prev[j-1],d));
            cur[j]=fmax(d,m);
        }
        { double *t=prev; prev=cur; cur=t; }
    }
    result=prev[nb-1];
    free(prev); free(cur);
    return result;
}
