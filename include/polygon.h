#ifndef GTESS_POLYGON_H
#define GTESS_POLYGON_H
#include <stddef.h>
#include "point.h"

int GT_Polygon_checkConvex(size_t n, const GT_Point points[static n], int cw);
int GT_Polygon_checkOrigin(size_t n, const GT_Point points[static n], int cw);
int GT_Polygon_contains(size_t n, const GT_Point points[static n], GT_Point p, int cw);
double GT_Polygon_area(size_t n, const GT_Point points[static n]);

double GT_Polygon_diameter(size_t n, const GT_Point points[static n], int cw);
size_t GT_Polygon_intersectConvex(GT_Point *points_c, size_t n_a, const GT_Point points_a[static n_a], size_t n_b, const GT_Point points_b[static n_b]);

#endif //GTESS_POLYGON_H

