#include <stddef.h>
#include <math.h>
#include "point.h"
#include "polygon.h"

int GT_Polygon_checkConvex(size_t n, const GT_Point points[static n], int cw){
	if(n < 3){
		return 1;
	}else for(size_t i = 0; i < n; ++i){
		if((cw ? -1 : 1)*GT_Point_cross_abc(points[(i + 1)%n], points[i], points[(n + i - 1)%n]) < GT_EPSILON){
			return 0;
		}
	}
	return 1;
}

int GT_Polygon_checkOrigin(size_t n, const GT_Point points[static n], int cw){
	if(n < 2){
		return !n || GT_Point_sqmag(points[0]) < GT_EPSILON;
	}else for(size_t i = 0; i < n; ++i){
		if((cw ? -1 : 1)*GT_Point_cross(points[i], points[(i + 1)%n]) < GT_EPSILON){
			return 0;
		}
	}
	return 1;
}

int GT_Polygon_contains(size_t n, const GT_Point points[static n], GT_Point p, int cw){
	if(n < 2){
		return !n || GT_Point_sqdist(points[0], p) < GT_EPSILON;
	}else for(size_t i = 0; i < n; ++i){
		if((cw ? -1 : 1)*GT_Point_cross_abc(points[i], p, points[(i + 1)%n]) < GT_EPSILON){
			return 0;
		}
	}
	return 1;
}

double GT_Polygon_area(size_t n, const GT_Point points[static n]){
    if(n < 3){
        return 0;
    }else if(n == 3){
        return GT_Point_cross_abc(points[2], points[1], points[0]) / 2;
    }
    double a = GT_Point_cross(points[n - 1], points[0]);
    for(size_t i = 1; i < n; ++i){
        a += GT_Point_cross(points[i - 1], points[i]);
    }
    return a/2;
}

double GT_Polygon_diameter(size_t n, const GT_Point points[static n], int cw){
	if(n < 3){
		switch(n){
		case 0:
			return INFINITY;
		case 1:
			return 0;
		case 2:
			return GT_Point_dist(points[0], points[1]);
		}
	}
	size_t top_i = 2;
	GT_Point bottom_axis = GT_Point_unit(GT_Point_sub(points[1], points[0]));
	GT_Point caliper_axis = GT_Point_ccw(bottom_axis);
	double height = GT_Point_dot(caliper_axis, points[top_i]);
	for(size_t i = 3; i < n; ++i){
		double height_here = GT_Point_dot(caliper_axis, points[i]);
		if(height_here < height){
			break;
		}
		height = height_here;
		top_i = i;
	}
	if(fabs(height) < GT_EPSILON){//the points are colinear
		double min_dot = GT_Point_dot(bottom_axis, points[0]);
		double max_dot = GT_Point_dot(bottom_axis, points[1]);
		for(size_t i = 2; i < n; ++i){
			double dot = GT_Point_dot(bottom_axis, points[i]);
			if(dot < max_dot){
				break;
			}
			max_dot = dot;
		}
		for(size_t i = n - 1; i > 1; --i){
			double dot = GT_Point_dot(bottom_axis, points[i]);
			if(dot > min_dot){
				break;
			}
			min_dot = dot;
		}
		return max_dot - min_dot;
	}
	double sqdiam = GT_Point_sqdist(points[0], points[top_i]);
	for(size_t bottom_i = 1; bottom_i < n; ++bottom_i){
		bottom_axis = GT_Point_unit(GT_Point_sub(points[(bottom_i + 1)%n], points[bottom_i]));
		caliper_axis = GT_Point_ccw(bottom_axis);
		height = GT_Point_dot(caliper_axis, points[top_i]);
		for(size_t i = (top_i + 1)%n; i != bottom_i; i = (i + 1)%n){
			double height_here = GT_Point_dot(caliper_axis, points[i]);
			if(height_here < height){
				break;
			}
			height = height_here;
			top_i = i;
		}
		double sqdiam_here = GT_Point_sqdist(points[bottom_i], points[top_i]);
		if(sqdiam_here > sqdiam){
			sqdiam = sqdiam_here;
		}
	}
	return sqrt(sqdiam);
}

