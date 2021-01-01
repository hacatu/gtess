#ifndef GTESS_POLYGON_H
#define GTESS_POLYGON_H
#include <stddef.h>
#include <inttypes.h>
#include "point.h"

typedef struct{
	GT_Point *points;
	size_t len, cap;
} GT_Polygon;

typedef struct{
	GT_Point point, in_normal;
} GT_Halfplane;

int GT_Polygon_checkConvex(size_t n, const GT_Point points[static n], int cw);
int GT_Polygon_checkOrigin(size_t n, const GT_Point points[static n], int cw);
int GT_Polygon_checkOrigin_strided(size_t n, size_t stride, const void *buf, int cw);
int GT_Polygon_contains(size_t n, const GT_Point points[static n], GT_Point p, int cw);
int GT_Halfplane_contains(GT_Halfplane halfplane, GT_Point p);
double GT_Polygon_area(size_t n, const GT_Point points[static n]);

size_t GT_Polygon_nearestVertexExternal(size_t n, const GT_Point points[static n], GT_Point a);

size_t GT_Polygon_properVertices(GT_Point *out, size_t n, const GT_Point points[static n]);

int GT_Polygon_equal_nubbed(size_t n_a, const GT_Point points_a[static n_a], size_t n_b, const GT_Point points_b[static n_b]);
int GT_Polygon_equal_stacked(size_t n_a, const GT_Point points_a[static n_a], size_t n_b, const GT_Point points_b[static n_b]);
int GT_Halfplane_equal(GT_Halfplane a, GT_Halfplane b);

double GT_Polygon_diameter(size_t n, const GT_Point points[static n], size_t *ai, size_t *bi, int cw);
size_t GT_Polygon_dimension(size_t n, const GT_Point points[static n]);
GT_Point GT_Polygon_interiorPoint(size_t n, const GT_Point points[static n]);
GT_Point GT_Polygon_centroid(size_t n, const GT_Point points[static n]);
size_t GT_Polygon_intersectConvex(GT_Point *points_c, size_t n_a, const GT_Point points_a[static n_a], size_t n_b, const GT_Point points_b[static n_b]);
size_t GT_Polygon_intersectConvexHalfplane(GT_Point *out, size_t n, const GT_Point points[static n], GT_Halfplane halfplane);
int GT_Halfplane_intersectHalfplane(GT_Point *out, GT_Halfplane a, GT_Halfplane b);
int GT_Halfplane_intersectSegment(GT_Point *out, GT_Halfplane halfplane, GT_Point a, GT_Point b);

uint64_t GT_gcd(uint64_t a, uint64_t b);

#endif //GTESS_POLYGON_H

