#define _GNU_SOURCE
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "point.h"
#include "polygon.h"
#include "polygon_generators.h"

#define DEBUG_PRINTF(do_print, printf_args...) (do_print ? printf(printf_args) : 0)

static int checkDiamondIntersections(int print_err){
	GT_Point points_a[] = {{1, 2}, {0, 3}, {-1, 2}, {0, 1}};
	GT_Point points_b[] = {{4, 2}, {3, 3}, {2, 2}, {3, 1}};
	GT_Point points_c[8];
	size_t n_c = GT_Polygon_intersectConvex(points_c, NULL, 4, points_a, 4, points_b);
	if(n_c){
		DEBUG_PRINTF(print_err, "Intersection found between horizontally separated polygons.\n");
		return 0;
	}
	memcpy(points_b, (GT_Point[4]){{-1, -2}, {0, -3}, {1, -2}, {0, -1}}, 4*sizeof(GT_Point));
	n_c = GT_Polygon_intersectConvex(points_c, NULL, 4, points_a, 4, points_b);
	if(n_c){
		DEBUG_PRINTF(print_err, "Intersection found between vertically separated polygons.\n");
		return 0;
	}
	memcpy(points_b, (GT_Point[4]){{1, 1}, {0, 2}, {-1, 1}, {0, 0}}, 4*sizeof(GT_Point));
	n_c = GT_Polygon_intersectConvex(points_c, NULL, 4, points_a, 4, points_b);
	if(n_c != 4){
		DEBUG_PRINTF(print_err, "Intersection between half vertically overlapping diamonds didn't have 4 points.\n");
		return 0;
	}
	if(!GT_Polygon_equal_stacked(n_c, points_c, 4, (GT_Point[4]){{.5, 1.5}, {0, 2}, {-.5, 1.5}, {0, 1}})){
		DEBUG_PRINTF(print_err, "Intersection between half vertically overlapping diamonds was not correct.\n");
		return 0;
	}
	return 1;
}

int main(){
	if(!checkDiamondIntersections(1)){
		exit(EXIT_FAILURE);
	}
	printf("Passed polygon intersection sanity checks.\n");
}

