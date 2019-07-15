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

int GT_Halfplane_contains(GT_Halfplane halfplane, GT_Point p){
	return GT_Point_dot(halfplane.in_normal, GT_Point_sub(p, halfplane.point)) >= 0;
}

int GT_Polygon_equal_nubbed(size_t n_a, const GT_Point points_a[static n_a], size_t n_b, const GT_Point points_b[static n_b]){
	if(!n_a){
		return !n_b;
	}else if(!n_b){
		return 0;
	}
	size_t k = 0;
	for(; k < n_b; ++k){
		if(GT_Point_sqdist(points_a[0], points_b[k]) < GT_EPSILON){
			break;
		}
	}
	if(k == n_b){
		return 0;
	}
	for(; n_a > 1; --n_a){
		if(GT_Point_sqdist(points_a[0], points_a[n_a - 1]) >= GT_EPSILON){
			break;
		}
	}
	while(1){
		size_t a_i = 0, a_j = 1, b_i = k, b_j = (k + 1)%n_b;
		for(; a_j < n_a; ++a_j){
			if(GT_Point_sqdist(points_a[a_i], points_a[a_j]) >= GT_EPSILON){
				break;
			}
		}
		a_i = a_j;
		for(; b_j != k; b_j = (b_j + 1)%n_b){
			if(GT_Point_sqdist(points_b[b_i], points_b[b_j]) >= GT_EPSILON){
				break;
			}
		}
		b_i = b_j;
		if(a_i == n_a){
			return b_i == k;
		}else if(b_i == k){
			return 0;
		}else if(GT_Point_sqdist(points_a[a_i], points_b[b_i]) >= GT_EPSILON){
			return 0;
		}
	}
}

int GT_Polygon_equal_stacked(size_t n_a, const GT_Point points_a[static n_a], size_t n_b, const GT_Point points_b[static n_b]){
	if(n_a != n_b){
		return 0;
	}else if(!n_a){
		return 1;
	}
	size_t a_i = 0;
	if(GT_Point_sqdist(points_a[0], points_a[n_a - 1]) < GT_EPSILON){
		for(; ++a_i < n_a;){
			if(GT_Point_sqdist(points_a[0], points_a[a_i]) >= GT_EPSILON){
				break;
			}
		}
		for(size_t i = 0; i < n_b; ++i){
			if(GT_Point_sqdist(points_a[0], points_b[i]) >= GT_EPSILON){
				return 0;
			}
		}
		return 1;
	}
	size_t k = 0;
	for(; k < n_b; ++k){
		if(GT_Point_sqdist(points_a[a_i], points_b[k]) < GT_EPSILON){
			break;
		}
	}
	if(k == n_b){
		return 0;
	}
	for(size_t i = 1; i < n_a; ++i){
		if(GT_Point_sqdist(points_a[(a_i + i)%n_a], points_b[(k + i)%n_b]) >= GT_EPSILON){
			return 0;
		}
	}
	return 1;
}

int GT_Halfplane_equal(GT_Halfplane a, GT_Halfplane b){
	if(GT_Point_dot(a.in_normal, b.in_normal) < 1 - GT_EPSILON){
		return 0;
	}
	return fabs(GT_Point_dot(a.in_normal, GT_Point_sub(a.point, b.point))) <= GT_EPSILON;
}

size_t GT_Polygon_nearestVertexExternal(size_t n, const GT_Point points[static n], GT_Point a){
	if(n < 2){
		return 0;
	}
	size_t min_i = 0;
	double min_sqdist = GT_Point_sqdist(points[0], a);
	double sqdist = GT_Point_sqdist(points[1], a);
	if(sqdist <= min_sqdist){
		min_sqdist = sqdist;
		min_i = 1;
		for(size_t i = 2; i < n; ++i){
			sqdist = GT_Point_sqdist(points[i], a);
			if(sqdist > min_sqdist){
				break;
			}
			min_sqdist = sqdist;
			min_i = i;
		}
	}else for(size_t i = n - 1; i; --i){
		sqdist = GT_Point_sqdist(points[i], a);
		if(sqdist > min_sqdist){
			break;
		}
		min_sqdist = sqdist;
		min_i = i;
	}
	return min_i;
}

size_t GT_Polygon_properVertices(GT_Point *out, size_t n, const GT_Point points[static n]){
	if(!n){
		return 0;
	}
	size_t ai = n - 1, ci = 1, bi = 0;
	while(GT_Point_sqdist(points[ai], points[bi]) < GT_EPSILON){
		if(ai-- == 1){
			out[0] = points[0];
			return 1;
		}
	}
	while(GT_Point_sqdist(points[ci], points[bi]) < GT_EPSILON){
		if(++ci == ai){
			out[0] = points[0];
			out[1] = points[ai];
			return 2;
		}
	}
	size_t final_i = ai;
	size_t j = 0;
	while(1){
		if(GT_Point_dot_abc(points[ai], points[bi], points[ci]) >= GT_EPSILON){
			out[j++] = points[bi];
		}
		if(ci == final_i){
			ai = bi;
			bi = ci;
			ci = 0;
		}else if(bi == final_i){
			break;
		}else{
			size_t ni = ci + 1;
			while(1){
				if(GT_Point_sqdist(points[ci], points[ni]) >= GT_EPSILON){
					ai = bi;
					bi = ci;
					ci = ni;
					break;
				}else if(ni == final_i){
					ai = bi;
					bi = ci;
					ci = 0;
					break;
				}
			}
		}
	}
	return j;
}

size_t GT_Polygon_dimension(size_t n, const GT_Point points[static n]){
	if(!n){
		return 0;
	}
	GT_Point unique[3] = {[0]=points[0]};
	size_t dim = 0, i = 1;
	for(; i < n; ++i){
		if(GT_Point_cmp_xy(points[i], unique[0])){
			unique[1] = points[i];
			dim = 1;
			break;
		}
	}
	for(; i < n; ++i){
		if(GT_Point_cmp_xy(points[i], unique[0]) && GT_Point_cmp_xy(points[i], unique[1])){
			if(GT_Point_cross_abc(points[i], unique[0], unique[1])){
				return 2;
			}
		}
	}
	return dim;
}

double GT_Polygon_area(size_t n, const GT_Point points[static n]){
    if(n < 3){
        return 0;
    }else if(n == 3){
        return GT_Point_cross_abc(points[2], points[1], points[0])/2;
    }
    double a = GT_Point_cross(points[n - 1], points[0]);
    for(size_t i = 1; i < n; ++i){
        a += GT_Point_cross(points[i - 1], points[i]);
    }
    return a/2;
}

double GT_Polygon_diameter(size_t n, const GT_Point points[static n], size_t *ai, size_t *bi, int cw){
	if(n < 3){
		switch(n){
		case 0:
			return INFINITY;
		case 1:
			return 0;
		case 2:
			if(ai){
				*ai = 0;
				*bi = 1;
			}
			return GT_Point_dist(points[0], points[1]);
		}
	}
	size_t top_i = 2;
	GT_Point bottom_axis = GT_Point_unit(GT_Point_sub(points[1], points[0]));
	GT_Point caliper_axis = (cw ? GT_Point_cw : GT_Point_ccw)(bottom_axis);
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
		if(ai){
			*ai = 0;
			*bi = 1;
		}
		for(size_t i = 2; i < n; ++i){
			double dot = GT_Point_dot(bottom_axis, points[i]);
			if(dot < max_dot){
				break;
			}
			max_dot = dot;
			if(ai){
				*bi = i;
			}
		}
		for(size_t i = n - 1; i > 1; --i){
			double dot = GT_Point_dot(bottom_axis, points[i]);
			if(dot > min_dot){
				break;
			}
			min_dot = dot;
			if(ai){
				*ai = i;
			}
		}
		return max_dot - min_dot;
	}
	double sqdiam = GT_Point_sqdist(points[0], points[top_i]);
	if(ai){
		*ai = 0;
		*bi = top_i;
	}
	for(size_t bottom_i = 1; bottom_i < n; ++bottom_i){
		bottom_axis = GT_Point_unit(GT_Point_sub(points[(bottom_i + 1)%n], points[bottom_i]));
		caliper_axis = (cw ? GT_Point_cw : GT_Point_ccw)(bottom_axis);
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
			if(ai){
				*ai = top_i;
				*bi = bottom_i;
			}
		}
	}
	return sqrt(sqdiam);
}

uint64_t GT_gcd(uint64_t a, uint64_t b){
	if(a < b){
		b %= a;
	}
	while(1){
		if(!b){
			return a;
		}
		a %= b;
		if(!a){
			return b;
		}
		b %= a;
	}
}

