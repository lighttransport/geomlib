/*
 * geomlib - Convex polygon clipping (see geomlib_polygon_clip.h)
 *
 * Clean-room Sutherland-Hodgman edge-by-edge clipping.
 * C11, no external dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polygon_clip.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double side(geomlib_vec2 a, geomlib_vec2 b, geomlib_vec2 p)
{
    return ((double)b.x-a.x)*((double)p.y-a.y)
         - ((double)b.y-a.y)*((double)p.x-a.x);
}

int geomlib_clip_polygon_convex(const geomlib_vec2 *subject, size_t ns,
                                const geomlib_vec2 *clip, size_t nc,
                                geomlib_polygon2d *out)
{
    geomlib_vec2 *cur=NULL,*next=NULL;
    size_t ncur=0,nnext=0;
    size_t i,j;
    int rc=GEOMLIB_MESH_ERROR_INVALID;

    if (!out || !subject || !clip || ns<3 || nc<3) return rc;
    memset(out,0,sizeof(*out));
    cur=(geomlib_vec2*)malloc((ns+nc+8)*sizeof(geomlib_vec2));
    next=(geomlib_vec2*)malloc((ns+nc+8)*sizeof(geomlib_vec2));
    if (!cur||!next) goto cleanup;
    memcpy(cur,subject,ns*sizeof(geomlib_vec2));
    ncur=ns;

    for (j=0;j<nc && ncur;++j) {
        geomlib_vec2 a=clip[j],b=clip[(j+1)%nc];
        nnext=0;
        for (i=0;i<ncur;++i) {
            geomlib_vec2 s=cur[i],e=cur[(i+1)%ncur];
            double ds=side(a,b,s),de=side(a,b,e);
            int ins=ds>=-1e-12, ine=de>=-1e-12;
            if (ins) next[nnext++]=s;
            if (ins!=ine) {
                double t=ds/(ds-de);
                next[nnext++]=geomlib_vec2_create(
                    (float)(s.x+t*(e.x-s.x)),(float)(s.y+t*(e.y-s.y)));
            }
        }
        { geomlib_vec2 *tmp=cur; cur=next; next=tmp; }
        ncur=nnext;
    }
    if (ncur>2) {
        if (geomlib_polygon2d_init(out,ncur)!=GEOMLIB_MESH_SUCCESS) goto cleanup;
        for (i=0;i<ncur;++i)
            if (geomlib_polygon2d_add_vertex(out,cur[i])!=GEOMLIB_MESH_SUCCESS)
                goto cleanup;
    }
    rc=GEOMLIB_MESH_SUCCESS;

cleanup:
    free(cur); free(next);
    if (rc!=GEOMLIB_MESH_SUCCESS) geomlib_polygon2d_free(out);
    return rc;
}
