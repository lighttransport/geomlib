/*
 * geomlib - OBB / OBB overlap
 *
 * Separating-axis test over two oriented boxes (15 axes).
 */

#ifndef GEOMLIB_OBB_OBB_H
#define GEOMLIB_OBB_OBB_H

#include "geomlib_collision.h"

#ifdef __cplusplus
extern "C" {
#endif

int geomlib_obb_obb_overlap(const geomlib_obb *a, const geomlib_obb *b);

#ifdef __cplusplus
}
#endif

#endif /* GEOMLIB_OBB_OBB_H */
