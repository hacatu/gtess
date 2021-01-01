#define _GNU_SOURCE
#include <stddef.h>
#include <stdarg.h>
#include <math.h>
#include "point.h"

//basic arithmetic functions
GT_Point GT_Point_add(GT_Point a, GT_Point b){
    return (GT_Point){a.x + b.x, a.y + b.y};
}

GT_Point GT_Point_sub(GT_Point a, GT_Point b){
    return (GT_Point){a.x - b.x, a.y - b.y};
}

GT_Point GT_Point_scale(GT_Point x, double a){
    return (GT_Point){a*x.x, a*x.y};
}

GT_Point GT_Point_neg(GT_Point a){
    return (GT_Point){-a.x, -a.y};
}

//dot product functions
double GT_Point_dot(GT_Point a, GT_Point b){
    return a.x*b.x + a.y*b.y;
}

double GT_Point_dot_abc(GT_Point a, GT_Point b, GT_Point c){
    return (a.x - b.x)*(c.x - b.x) + (a.y - b.y)*(c.y - b.y);
}

double GT_Point_cosangle(GT_Point a){
    return a.x/GT_Point_mag(a);
}

double GT_Point_cosangle_between(GT_Point a, GT_Point b){
    return GT_Point_dot(a, b)/(GT_Point_mag(a)*GT_Point_mag(b));
}

double GT_Point_cosangle_abc(GT_Point a, GT_Point b, GT_Point c){
    return GT_Point_dot_abc(a, b, c)/(GT_Point_dist(a, b)*GT_Point_dist(c, b));
}

//cross product functions
double GT_Point_cross(GT_Point a, GT_Point b){
    return a.x*b.y - a.y*b.x;
}

double GT_Point_cross_abc(GT_Point a, GT_Point b, GT_Point c){
    return (a.x - b.x)*(c.y - b.y) - (a.y - b.y)*(c.x - b.x);
}

double GT_Point_sinangle(GT_Point a){
    return a.y/GT_Point_mag(a);
}

double GT_Point_sinangle_between(GT_Point a, GT_Point b){
    return GT_Point_cross(a, b)/(GT_Point_mag(a)*GT_Point_mag(b));
}

double GT_Point_sinangle_abc(GT_Point a, GT_Point b, GT_Point c){
    return GT_Point_cross_abc(a, b, c)/(GT_Point_dist(a, b)*GT_Point_dist(c, b));
}

double GT_Point_area(size_t n, ...){
    if(n < 3){
        return 0;
    }
    va_list args;
    va_start(args, n);
    if(n == 3){
        GT_Point a = va_arg(args, GT_Point), b = va_arg(args, GT_Point), c = va_arg(args, GT_Point);
        va_end(args);
        return GT_Point_cross_abc(c, b, a)/2;
    }
    GT_Point a0 = va_arg(args, GT_Point), a, b = a0;
    double A = 0;
    for(size_t i = 1; i < n; ++i){
        a = b;
        b = va_arg(args, GT_Point);
        A += GT_Point_cross(a, b);
    }
    A += GT_Point_cross(b, a0);
    return A/2;
}

//angle functions
double GT_Point_angle(GT_Point a){
    const double ret = atan2(a.y, a.x);
    return ret < 0 ? ret + 2*M_PI : ret;
}

double GT_Point_angle_between(GT_Point a, GT_Point b){
    const double ret = atan2(b.y, b.x) - atan2(a.y, a.x);
    return ret < 0 ? ret + 2*M_PI : ret;
}

double GT_Point_angle_abc(GT_Point a, GT_Point b, GT_Point c){
    const double ret = atan2(c.y - b.y, c.x - b.x) - atan2(a.y - b.y, a.x - b.x);
    return ret < 0 ? ret + 2*M_PI : ret;
}

//distance and magnitude functions
double GT_Point_dist(GT_Point a, GT_Point b){
    return hypot(b.x - a.x, b.y - a.y);
}

double GT_Point_sqdist(GT_Point a, GT_Point b){
    const double dx = b.x - a.x, dy = b.y - a.y;
    return dx*dx + dy*dy;
}

double GT_Point_mag(GT_Point a){
    return hypot(a.x, a.y);
}

double GT_Point_sqmag(GT_Point a){
    return a.x*a.x + a.y*a.y;
}

GT_Point GT_Point_unit(GT_Point a){
    const double mag = GT_Point_mag(a);
    return (GT_Point){a.x/mag, a.y/mag};
}

//rotation functions
GT_Point GT_Point_cw(GT_Point a) {
    return (GT_Point){a.y, -a.x};
}

GT_Point GT_Point_ccw(GT_Point a){
    return (GT_Point){-a.y, a.x};
}

GT_Point GT_Point_rotate(GT_Point x, double a){
    const double ca = cos(a), sa = sin(a);
    return (GT_Point){ca*x.x - sa*x.y, sa*x.x + ca*x.y};
}

GT_Point GT_Point_mirror(GT_Point x, GT_Point m){
    const double dx = x.x - m.x, dy = x.y - m.y;
    return (GT_Point){m.x - dx, m.y - dy};
}

GT_Point GT_Point_cw_around(GT_Point a, GT_Point o){
    const double dx = a.x - o.x, dy = a.y - o.y;
    return (GT_Point){o.x + dy, o.y - dx};
}

GT_Point GT_Point_ccw_around(GT_Point a, GT_Point o){
    const double dx = a.x - o.x, dy = a.y - o.y;
    return (GT_Point){o.x - dy, o.y + dx};
}

GT_Point GT_Point_rotate_around(GT_Point x, GT_Point o, double a){
    const double dx = x.x - o.x, dy = x.y - o.y;
    const double ca = cos(a), sa = sin(a);
    return (GT_Point){o.x + ca*dx - sa*dy, o.y + sa*dx + ca*dy};
}

//linear combination functions
GT_Point GT_Point_midpoint(GT_Point a, GT_Point b){
    return (GT_Point){(a.x + b.x)/2, (a.y + b.y)/2};
}

GT_Point GT_Point_centroid(size_t n, ...){
    va_list args;
    va_start(args, n);
    double x = 0, y = 0;
    for(size_t i = 0; i < n; ++ i){
        GT_Point a = va_arg(args, GT_Point);
        x += a.x;
        y += a.y;
    }
    va_end(args);
    return (GT_Point){x/n, y/n};
}

GT_Point GT_Point_lincomb2(GT_Point x, double a, GT_Point y, double b){
    return (GT_Point){x.x*a + y.x*b, x.y*a + y.y*b};
}

GT_Point GT_Point_lincomb(size_t n, ...){
    va_list args;
    va_start(args, n);
    double x = 0, y = 0;
    for(size_t i = 0; i < n; ++i){
        GT_Point a = va_arg(args, GT_Point);
        double w = va_arg(args, double);
        x += a.x*w;
        y += a.y*w;
    }
    va_end(args);
    return (GT_Point){x, y};
}

//projection functions
GT_Point GT_Point_project(GT_Point a, GT_Point b){
    return GT_Point_scale(b, GT_Point_dot(a, b)/GT_Point_sqmag(b));
}

GT_Point GT_Point_project_unit(GT_Point a, GT_Point b){
    return GT_Point_scale(b, GT_Point_dot(a, b));
}

GT_Point GT_Point_project_sqmag(GT_Point a, GT_Point b, double b_sqmag){
    return GT_Point_scale(b, GT_Point_dot(a, b)/b_sqmag);
}

//line intersection functions
int GT_Point_intersect_segments_open_ab(GT_Point *out, GT_Point a_a, GT_Point a_b, GT_Point b_a, GT_Point b_b){
	return GT_Point_intersect_segments_open_po(out, a_a, GT_Point_sub(a_b, a_a), b_a, GT_Point_sub(b_b, b_a));
}

int GT_Point_intersect_segments_open_po(GT_Point *out, GT_Point a_p, GT_Point a_o, GT_Point b_p, GT_Point b_o){
	double det = GT_Point_cross(a_o, b_o);
	if(fabs(det) < GT_EPSILON){
		if(fabs(GT_Point_cross(GT_Point_sub(b_p, a_p), a_o)) < GT_EPSILON){//collinear
			double rsqmag = GT_Point_sqmag(a_o);
			double t0 = GT_Point_dot(GT_Point_sub(b_p, a_p), a_o)/rsqmag;
			double t1 = t0 + GT_Point_dot(b_o, a_o)/rsqmag;
			if(t1 < t0){
				if(t1 <= 0){
					if(t0 > 0){
						*out = a_p;
						return 2;
					}
				}else if(t1 < 1){
					*out = GT_Point_add(a_p, a_o);
					return 2;
				}
			}else if(t0 <= 0){
				if(t1 > 0){
					*out = a_p;
					return 2;
				}
			}else if(t0 < 1){
				*out = b_p;
				return 2;
			}
		}
		return 0;
	}
	double t = GT_Point_cross(GT_Point_sub(b_p, a_p), a_o)/det;
	double u = GT_Point_cross(GT_Point_sub(b_p, a_p), b_o)/det;
	if(0 < t && t < 1 && 0 < u && u < 1){
		*out = GT_Point_add(a_p, GT_Point_scale(a_o, t));
		return 1;
	}
	return 0;
}

int GT_Point_intersect_segments_closed_ab(GT_Point *out, GT_Point a_a, GT_Point a_b, GT_Point b_a, GT_Point b_b){
	return GT_Point_intersect_segments_closed_po(out, a_a, GT_Point_sub(a_b, a_a), b_a, GT_Point_sub(b_b, b_a));
}

int GT_Point_intersect_segments_closed_po(GT_Point *out, GT_Point a_p, GT_Point a_o, GT_Point b_p, GT_Point b_o){
	double det = GT_Point_cross(a_o, b_o);
	if(fabs(det) < GT_EPSILON){
		if(fabs(GT_Point_cross(GT_Point_sub(b_p, a_p), a_o)) < GT_EPSILON){//collinear
			double rsqmag = GT_Point_sqmag(a_o);
			double t0 = GT_Point_dot(GT_Point_sub(b_p, a_p), a_o)/rsqmag;
			double t1 = t0 + GT_Point_dot(b_o, a_o)/rsqmag;
			if(t1 < t0){
				if(t1 <= 0){
					if(t0 >= 0){
						*out = a_p;
						return 2;
					}
				}else if(t1 <= 1){
					*out = GT_Point_add(a_p, a_o);
					return 2;
				}
			}else if(t0 <= 0){
				if(t1 >= 0){
					*out = a_p;
					return 2;
				}
			}else if(t0 <= 1){
				*out = b_p;
				return 2;
			}
		}
		return 0;
	}
	double t = GT_Point_cross(GT_Point_sub(b_p, a_p), a_o)/det;
	double u = GT_Point_cross(GT_Point_sub(b_p, a_p), b_o)/det;
	if(0 <= t && t <= 1 && 0 <= u && u <= 1){
		*out = GT_Point_add(a_p, GT_Point_scale(a_o, t));
		return 1;
	}
	return 0;
}

double GT_Point_slope(GT_Point a, GT_Point b){
	return (b.x - a.x)/(b.y - a.y);
}

int GT_Point_cmp_xy(GT_Point a, GT_Point b){
	if(a.x < b.x){
		return -1;
	}else if(a.x > b.x){
		return 1;
	}else if(a.y < b.y){
		return -1;
	}else if(a.y > b.y){
		return 1;
	}
	return 0;
}

int GT_Point_cmp_x(GT_Point a, GT_Point b){
	if(a.x < b.x){
		return -1;
	}else if(a.x > b.x){
		return 1;
	}
	return 0;
}

int GT_Point_cmp_fn_x(const void *a, const void *b){
	return GT_Point_cmp_x(*(const GT_Point*)a, *(const GT_Point*)b);
}

GT_Point GT_Point_e1 = {1, 0}, GT_Point_e2 = {0, 1}, GT_Point_zero = {0, 0}, GT_Point_neg1_polar = {M_PI, 1};

GT_Point GT_Point_dir[4] = {
	[GT_Direction_right]={1, 0},
	[GT_Direction_up]={0, 1},
	[GT_Direction_left]={-1, 0},
	[GT_Direction_down]={0, -1}
};

