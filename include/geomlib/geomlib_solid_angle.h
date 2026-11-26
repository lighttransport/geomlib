/*
 * geomlib - Solid angle and winding number
 *
 * Exact signed solid angle subtended by triangles, polygons and closed meshes
 * at a query point, using the Van Oosterom-Strackee formula.
 */

#ifndef GEOMLIB_SOLID_ANGLE_H
#define GEOMLIB_SOLID_ANGLE_H

#include "geomlib_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Signed solid angle (steradians) of triangle abc seen from p, in [-2pi,2pi].
 * The sign follows the winding of (a,b,c) as viewed from p. */
double geomlib_solid_angle_triangle(geomlib_vec3 p, geomlib_vec3 a,
                                    geomlib_vec3 b, geomlib_vec3 c);

/* Signed solid angle of a polygon seen from p (convex polygons, or any simple
 * polygon when p is coplanar-safe; uses a triangle fan from vertex 0). */
double geomlib_solid_angle_polygon(const geomlib_vec3 *verts, size_t n,
                                   geomlib_vec3 p);

/* Total signed solid angle of a closed mesh at p: 4pi inside, 0 outside for a
 * consistently outward-oriented mesh. */
double geomlib_solid_angle_mesh(const geomlib_trimesh *mesh, geomlib_vec3 p);

/* Winding number of a closed mesh at p, i.e. solid_angle / (4 pi). */
double geomlib_winding_number_mesh(const geomlib_trimesh *mesh, geomlib_vec3 p);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_SOLID_ANGLE_H */
