#define _GNU_SOURCE
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "point.h"

#define NUM_TRIALS 100

static void randomizePoints(gsl_rng *rng, size_t n, GT_Point points[static n]){
    for(size_t i = 0; i < n; ++i){
        points[i].x = gsl_ran_gaussian(rng, 100);
        points[i].y = gsl_ran_gaussian(rng, 100);
    }
}

int main(){
    gsl_rng_env_setup();
    const gsl_rng_type *rng_T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(rng_T);
    size_t passed_tests = 0;
    for(size_t i = 0; i < NUM_TRIALS; ++i){
        GT_Point p;
        randomizePoints(rng, 1, &p);
        double angle = GT_Point_angle(p);
        double fwd_angle = GT_Point_angle_between(GT_Point_e1, p);
        double rev_angle = GT_Point_angle_between(p, GT_Point_e1);
        double h = GT_Point_mag(p);
        GT_Point cw = GT_Point_cw(p);
        GT_Point u = GT_Point_unit(p);
        if(fabs(angle - fwd_angle) >= GT_EPSILON){
            continue;
        }else if(fabs(fwd_angle + rev_angle - 2*M_PI) >= GT_EPSILON){
            continue;
        }else if(fabs(GT_Point_sqmag(p) - GT_Point_sqdist(p, GT_Point_zero)) >= GT_EPSILON){
            continue;
        }else if(GT_Point_sqdist(cw, GT_Point_rotate(p, -M_PI/2)) >= GT_EPSILON){
            continue;
        }else if(GT_Point_sqdist(GT_Point_neg(cw), GT_Point_ccw(p)) >= GT_EPSILON){
            continue;
        }else if(fabs(p.y/h - GT_Point_sinangle(p)) >= GT_EPSILON){
            continue;
        }else if(fabs(p.x/h - GT_Point_cosangle(p)) >= GT_EPSILON){
            continue;
        }else if(fabs(angle - GT_Point_angle(u)) >= GT_EPSILON || fabs(GT_Point_mag(u) - 1) >= GT_EPSILON){
            continue;
        }
        ++passed_tests;
    }
    printf("Passed %zu/%zu 1 point tests\n", passed_tests, (size_t)(NUM_TRIALS));
    passed_tests = 0;
    for(size_t i = 0; i < NUM_TRIALS; ++i){
        GT_Point points[2];
        randomizePoints(rng, 2, points);
        double fwd_angle = GT_Point_angle_between(points[0], points[1]);
        double rev_angle = GT_Point_angle_between(points[1], points[0]);
        GT_Point cw = GT_Point_cw_around(points[0], points[1]);
        double a = GT_Point_mag(points[0]), b = GT_Point_mag(points[1]), c = GT_Point_dist(points[0], points[1]);
        double cos_a = GT_Point_cosangle_abc(GT_Point_zero, points[1], points[0]);
        double cos_b = GT_Point_cosangle_abc(points[1], points[0], GT_Point_zero);
        double cos_c = GT_Point_cosangle_between(points[0], points[1]);
        if(fabs(fwd_angle + rev_angle - 2*M_PI) >= GT_EPSILON) {
            continue;
        }else if(GT_Point_sqdist(cw, GT_Point_rotate_around(points[0], points[1], -M_PI/2)) >= GT_EPSILON){
            continue;
        }else if(GT_Point_sqdist(GT_Point_mirror(cw, points[1]), GT_Point_ccw_around(points[0], points[1])) >= GT_EPSILON){
            continue;
        }else if(fabs(a*a + b*b - c*c - 2*a*b*cos_c) >= GT_EPSILON || fabs(b*b + c*c - a*a - 2*b*c*cos_a) >= GT_EPSILON
                 || fabs(c*c + a*a - b*b - 2*c*a*cos_b) >= GT_EPSILON){
            continue;
        }else if(GT_Point_area(2, points[0], points[1])){
            continue;
        }
        ++passed_tests;
    }
    printf("Passed %zu/%zu 2 point tests\n", passed_tests, (size_t)(NUM_TRIALS));
    passed_tests = 0;
    for(size_t i = 0; i < NUM_TRIALS; ++i){
        GT_Point points[3];
        randomizePoints(rng, 3, points);
        double fwd_angle = GT_Point_angle_abc(points[0], points[1], points[2]);
        double rev_angle = GT_Point_angle_abc(points[2], points[1], points[0]);
        double b_angle = GT_Point_angle_abc(points[1], points[2], points[0]);
        double c_angle = GT_Point_angle_abc(points[2], points[0], points[1]);
        double interior = fwd_angle + b_angle + c_angle;
        if(interior >= M_PI + GT_EPSILON){
            interior = 6*M_PI - interior;
        }
        if(fabs(fwd_angle + rev_angle - 2*M_PI) >= GT_EPSILON){
            continue;
        }else if(fabs(interior - M_PI) >= GT_EPSILON){
            continue;
        }
        ++passed_tests;
    }
    printf("Passed %zu/%zu 3 point tests\n", passed_tests, (size_t)(NUM_TRIALS));
    gsl_rng_free(rng);
}

