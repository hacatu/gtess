#define _GNU_SOURCE
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "point.h"

#define NUM_TRIALS 100

static double initTriangle(GT_Point points[static 3], GT_Point o, double angle, double base_a, double base_b,
                           double directed_height){
    const double base = base_a + base_b;
    points[0] = o;
    const GT_Point point_opposite = GT_Point_add(o, GT_Point_rotate((GT_Point) {base, 0}, angle));
    const GT_Point point_apex = GT_Point_add(o, GT_Point_rotate((GT_Point) {base_a, directed_height}, angle));
    if(directed_height > 0){
        points[1] = point_opposite;
        points[2] = point_apex;
        return base*directed_height/2;
    }
    points[1] = point_apex;
    points[2] = point_opposite;
    return base*-directed_height/2;
}

static double randomTriangle(GT_Point points[static 3], gsl_rng *rng){
    const double x = gsl_ran_gaussian(rng, 100), y = gsl_ran_gaussian(rng, 100);
    const double angle = gsl_ran_flat(rng, 0, 2*M_PI);
    const double base_a = gsl_ran_exponential(rng, 5), base_b = gsl_ran_exponential(rng, 5);
    const double directed_height = gsl_ran_gaussian(rng, 10);
    return initTriangle(points, (GT_Point){x, y}, angle, base_a, base_b, directed_height);
}

int main(){
    gsl_rng_env_setup();
    const gsl_rng_type *rng_T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(rng_T);
    size_t passed_tests = 0;
    for(size_t i = 0; i < NUM_TRIALS; ++i){
        GT_Point points[3];
        randomTriangle(points, rng);
        GT_Point c1 = GT_Point_centroid(3, points[0], points[1], points[2]);
        GT_Point c2 = GT_Point_lincomb(3, points[0], 1./3, points[1], 1./3, points[2], 1./3);
        GT_Point a = GT_Point_sub(points[1], points[0]), b = GT_Point_sub(points[2], points[0]);
        GT_Point mid_ab = GT_Point_midpoint(a, b);
        double a_sqmag = GT_Point_sqmag(a);
        GT_Point b_proj_a = GT_Point_project_sqmag(b, a, a_sqmag);
        if(GT_Point_sqdist(c1, c2) >= GT_EPSILON){
            continue;
        }else if(GT_Point_sqdist(GT_Point_project_sqmag(mid_ab, a, a_sqmag), GT_Point_midpoint(b_proj_a, a)) >= GT_EPSILON){
            continue;
        }
        ++passed_tests;
    }
    printf("Passed %zu/%zu tests\n", passed_tests, (size_t)(NUM_TRIALS));
    gsl_rng_free(rng);
}

