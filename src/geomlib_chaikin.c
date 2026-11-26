/*
 * geomlib - Chaikin subdivision (see geomlib_chaikin.h)
 *
 * Clean-room corner cutting. C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_chaikin.h"

#include <stdlib.h>
#include <string.h>

size_t geomlib_chaikin_subdivide(const geomlib_vec3 *points, size_t count,
                                 int iterations, geomlib_vec3 *out,
                                 size_t out_capacity)
{
    geomlib_vec3 *cur=NULL,*next=NULL;
    size_t n=count,i,written=0;

    if (!points || count<2 || !out || !out_capacity || iterations<0) return 0;
    if ((count-1)*(1u<<(iterations>24?24:iterations))+1 > out_capacity)
        return 0;
    cur=(geomlib_vec3*)malloc(count*sizeof(geomlib_vec3));
    if (!cur) return 0;
    memcpy(cur,points,count*sizeof(geomlib_vec3));
    n=count;

    for (int it=0;it<iterations;++it) {
        size_t need=(n-1)*2+1;
        next=(geomlib_vec3*)realloc(next,need*sizeof(geomlib_vec3));
        if (!next) { free(cur); return 0; }
        for (i=0;i+1<n;++i) {
            next[2*i]=geomlib_vec3_add(geomlib_vec3_mul(cur[i],0.75f),
                                       geomlib_vec3_mul(cur[i+1],0.25f));
            next[2*i+1]=geomlib_vec3_add(geomlib_vec3_mul(cur[i],0.25f),
                                         geomlib_vec3_mul(cur[i+1],0.75f));
        }
        next[0]=cur[0];                 /* preserve open-curve endpoints */
        next[need-1]=cur[n-1];
        { geomlib_vec3 *t=cur; cur=next; next=t; }
        n=need;
    }
    for (i=0;i<n;++i) out[i]=cur[i];
    written=n;
    free(cur); free(next);
    return written;
}
