#ifndef GTESS_POINT_H
#define GTESS_POINT_H
#include <stddef.h>

#define GT_EPSILON 5e-6

typedef struct{
    double x, y;
} GT_Point;

typedef enum{
	GT_Direction_right,
	GT_Direction_up,
	GT_Direction_left,
	GT_Direction_down
} GT_Direction;

//basic arithmetic functions
GT_Point GT_Point_add(GT_Point a, GT_Point b);
GT_Point GT_Point_sub(GT_Point a, GT_Point b);
GT_Point GT_Point_scale(GT_Point x, double a);
GT_Point GT_Point_neg(GT_Point a);

//dot product functions
double GT_Point_dot(GT_Point a, GT_Point b);
double GT_Point_dot_abc(GT_Point a, GT_Point b, GT_Point c);
double GT_Point_cosangle(GT_Point a);
double GT_Point_cosangle_between(GT_Point a, GT_Point b);
double GT_Point_cosangle_abc(GT_Point a, GT_Point b, GT_Point c);

//cross product functions
double GT_Point_cross(GT_Point a, GT_Point b);
double GT_Point_cross_abc(GT_Point a, GT_Point b, GT_Point c);
double GT_Point_sinangle(GT_Point a);
double GT_Point_sinangle_between(GT_Point a, GT_Point b);
double GT_Point_sinangle_abc(GT_Point a, GT_Point b, GT_Point c);
double GT_Point_area(size_t n, ...);

//angle functions
double GT_Point_angle(GT_Point a);
double GT_Point_angle_between(GT_Point a, GT_Point b);
double GT_Point_angle_abc(GT_Point a, GT_Point b, GT_Point c);

//distance and magnitude functions
double GT_Point_dist(GT_Point a, GT_Point b);
double GT_Point_sqdist(GT_Point a, GT_Point b);
double GT_Point_mag(GT_Point a);
double GT_Point_sqmag(GT_Point a);
GT_Point GT_Point_unit(GT_Point a);

//rotation functions
GT_Point GT_Point_cw(GT_Point a);
GT_Point GT_Point_ccw(GT_Point a);
GT_Point GT_Point_rotate(GT_Point x, double a);
GT_Point GT_Point_mirror(GT_Point x, GT_Point m);
GT_Point GT_Point_cw_around(GT_Point a, GT_Point o);
GT_Point GT_Point_ccw_around(GT_Point a, GT_Point o);
GT_Point GT_Point_rotate_around(GT_Point x, GT_Point o, double a);

//linear combination functions
GT_Point GT_Point_midpoint(GT_Point a, GT_Point b);
GT_Point GT_Point_centroid(size_t n, ...);
GT_Point GT_Point_lincomb2(GT_Point x, double a, GT_Point y, double b);
GT_Point GT_Point_lincomb(size_t n, ...);

//projection functions
GT_Point GT_Point_project(GT_Point a, GT_Point b);
GT_Point GT_Point_project_unit(GT_Point a, GT_Point b);
GT_Point GT_Point_project_sqmag(GT_Point a, GT_Point b, double b_sqmag);

//line intersection functions
//endpoint a and endpoint b
int GT_Point_intersect_segments_open_ab(GT_Point *out, GT_Point a_a, GT_Point a_b, GT_Point b_a, GT_Point b_b);
//point p and offset o
int GT_Point_intersect_segments_open_po(GT_Point *out, GT_Point a_p, GT_Point a_o, GT_Point b_p, GT_Point b_o);
int GT_Point_intersect_segments_closed_ab(GT_Point *out, GT_Point a_a, GT_Point a_b, GT_Point b_a, GT_Point b_b);
int GT_Point_intersect_segments_closed_po(GT_Point *out, GT_Point a_p, GT_Point a_o, GT_Point b_p, GT_Point b_o);

double GT_Point_slope(GT_Point a, GT_Point b);

int GT_Point_cmp_xy(GT_Point a, GT_Point b);
int GT_Point_cmp_x(GT_Point a, GT_Point b);

int GT_Point_cmp_fn_x(const void *a, const void *b);

//common points
extern GT_Point GT_Point_e1, GT_Point_e2, GT_Point_zero, GT_Point_neg1_polar;
extern GT_Point GT_Point_dir[4];

#endif //GTESS_POINT_H

