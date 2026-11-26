/*
 * geomlib - Polygon moments (see geomlib_polygon_moments.h)
 *
 * Clean-room shoelace second-moment integration. C11, no external
 * dependencies, wasm32 safe.
 */

#include "geomlib/geomlib_polygon_moments.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int geomlib_compute_polygon_moments(const geomlib_vec2 *polygon, size_t n,
                            geomlib_polygon_moments *out)
{
    size_t i;
    double a2=0.0, cx=0.0, cy=0.0;
    double xx=0.0, xy=0.0, yy=0.0;
    double area, sign;

    if (!out || !polygon || n<3) return GEOMLIB_MESH_ERROR_INVALID;
    for (i=0;i<n;++i) {
        size_t j=(i+1)%n;
        double xi=polygon[i].x, yi=polygon[i].y;
        double xj=polygon[j].x, yj=polygon[j].y;
        double cr=xi*yj-xj*yi;
        a2+=cr;
        cx+=(xi+xj)*cr;
        cy+=(yi+yj)*cr;
        xx+=(yi*yi+yi*yj+yj*yj)*cr;
        yy+=(xi*xi+xi*xj+xj*xj)*cr;
        xy+=(2.0*xi*yi+xi*yj+xj*yi+2.0*xj*yj)*cr;
    }
    a2*=0.5;
    if (!isfinite(a2) || fabs(a2)<1e-18) return GEOMLIB_MESH_ERROR_INVALID;
    sign = a2>0.0 ? 1.0 : -1.0;
    area = fabs(a2);
    cx /= (6.0*a2);
    cy /= (6.0*a2);
    xx *= sign/12.0;
    yy *= sign/12.0;
    xy *= sign/24.0;
    /* Parallel-axis shift to centroid. */
    out->ixx=xx-area*cy*cy;
    out->iyy=yy-area*cx*cx;
    out->ixy=xy-area*cx*cy;
    out->area=area;
    out->centroid=geomlib_vec2_create((float)cx,(float)cy);
    return GEOMLIB_MESH_SUCCESS;
}
