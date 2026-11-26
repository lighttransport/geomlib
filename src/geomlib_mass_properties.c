/*
 * geomlib - Mass properties (see geomlib_mass_properties.h)
 *
 * Clean-room signed-tetrahedron integration. C11, no external dependencies,
 * wasm32 safe.
 */

#include "geomlib/geomlib_mass_properties.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_compute_mass_properties(const geomlib_trimesh *mesh,
                                    geomlib_mass_properties *out)
{
    size_t t, i, j;
    double vol=0.0;
    double com[3]={0,0,0};
    double second[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    double tr, d2;
    geomlib_vec3 c;

    if (!out || !mesh || !mesh->vertices || !mesh->indices ||
        !mesh->triangle_count) return GEOMLIB_MESH_ERROR_INVALID;
    memset(out,0,sizeof(*out));

    for (t=0;t<mesh->triangle_count;++t) {
        uint32_t ix[3]={mesh->indices[t*3],mesh->indices[t*3+1],mesh->indices[t*3+2]};
        geomlib_vec3 p[3];
        double v;
        double s[3]={0,0,0};
        for (i=0;i<3;++i) {
            if (ix[i]>=mesh->vertex_count) return GEOMLIB_MESH_ERROR_INVALID;
            p[i]=mesh->vertices[ix[i]];
            s[0]+=p[i].x; s[1]+=p[i].y; s[2]+=p[i].z;
        }
        v=geomlib_vec3_dot(p[0],geomlib_vec3_cross(p[1],p[2]))/6.0;
        for (i=0;i<3;++i) com[i]+=v*s[i]/4.0;
        /* Integral of x x^T over the tetrahedron (0,p0,p1,p2). */
        {
            {
                double q[3][3];
                for (i=0;i<3;++i) {
                    q[0][i]=(i==0)?p[0].x:((i==1)?p[1].x:p[2].x);
                    q[1][i]=(i==0)?p[0].y:((i==1)?p[1].y:p[2].y);
                    q[2][i]=(i==0)?p[0].z:((i==1)?p[1].z:p[2].z);
                }
                for (i=0;i<3;++i)
                    for (j=0;j<3;++j)
                        second[i][j]+=v/20.0*(q[i][0]*q[j][0]+q[i][1]*q[j][1]
                                             +q[i][2]*q[j][2]+s[i]*s[j]);
            }
        }
        vol+=v;
    }
    /* Tolerate consistently inverted winding by flipping all integrals. */
    if (vol<0.0) {
        vol=-vol;
        for (i=0;i<3;++i) com[i]=-com[i];
        for (i=0;i<3;++i) for (j=0;j<3;++j) second[i][j]=-second[i][j];
    }
    if (!(vol>1e-12) || !isfinite(vol)) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<3;++i) com[i]/=vol;
    c=geomlib_vec3_create((float)com[0],(float)com[1],(float)com[2]);

    tr=second[0][0]+second[1][1]+second[2][2];
    d2=com[0]*com[0]+com[1]*com[1]+com[2]*com[2];
    /* Inertia at origin = tr(I) - second; shift to the centre of mass. */
    for (i=0;i<3;++i)
        for (j=0;j<3;++j) {
            double v=((i==j)?tr:0.0)-second[i][j];
            v-=vol*(((i==j)?d2:0.0)-com[i]*com[j]);
            out->inertia.m[j*3+i]=(float)v;
        }
    out->volume=(float)vol;
    out->center=c;
    return GEOMLIB_MESH_SUCCESS;
}
