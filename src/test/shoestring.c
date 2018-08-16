#define _GNU_SOURCE
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "point.h"

#define NUM_TRIALS 100

static double computeAreaShoestring(size_t n, GT_Point points[static n]){
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

static double computeAreaTriangleSine(GT_Point points[static 3]){
    return GT_Point_sinangle_abc(points[1], points[0], points[2])
           *GT_Point_dist(points[1], points[0])*GT_Point_dist(points[2], points[0])/2;
}

static double computeAreaTriangleSine2(GT_Point points[static 3]){
    const GT_Point a = GT_Point_sub(points[1], points[0]);
    const GT_Point b = GT_Point_sub(points[2], points[0]);
    return GT_Point_sinangle_between(a, b)*GT_Point_mag(a)*GT_Point_mag(b)/2;
}

static double computeAreaTriangleProject(GT_Point points[static 3]){
    const GT_Point a = GT_Point_sub(points[1], points[0]);
    const GT_Point b = GT_Point_sub(points[2], points[1]);
    const GT_Point height = GT_Point_project(b, GT_Point_ccw(a));
    return GT_Point_cross(a, height)/2;
}

static int checkInteriorAngles(size_t n, GT_Point points[static n]){
    //the interior angles must sum to (n - 2)*pi
    double angle = 0;
    for(size_t i = 0; i < n; ++i){
        angle += GT_Point_angle_abc(points[(i + 1)%n], points[i], points[(i + n - 1)%n]);
    }
    //printf("%f | %f\n", (n - 2)*M_PI, angle);
    return fabs(angle - (n - 2)*M_PI) < GT_EPSILON;
}

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

static double initQuadrilateral(GT_Point points[static 4], GT_Point o, double angle, double base_a, double base_b,
                                double base_c, double directed_height_a, double directed_height_b){
    const double base = base_a + base_b + base_c;
    points[0] = o;
    const GT_Point point_opposite = GT_Point_add(o, GT_Point_rotate((GT_Point){base, 0}, angle));
    const GT_Point point_a = GT_Point_add(o, GT_Point_rotate((GT_Point){base_a, directed_height_a}, angle));
    const GT_Point point_b = GT_Point_add(o, GT_Point_rotate((GT_Point){base_a + base_b, directed_height_b}, angle));
    if(directed_height_a < 0){
        points[1] = point_a;
        if(directed_height_b < 0){
            points[2] = point_b;
            points[3] = point_opposite;
            return ((base_a + base_b)*-directed_height_a + (base_b + base_c)*-directed_height_b)/2;
        }
        points[2] = point_opposite;
        points[3] = point_b;
        return base*(-directed_height_a + directed_height_b)/2;
    }
    points[3] = point_a;
    if(directed_height_b < 0){
        points[1] = point_b;
        points[2] = point_opposite;
        return base*(directed_height_a + -directed_height_b)/2;
    }
    points[1] = point_opposite;
    points[2] = point_b;
    return ((base_a + base_b)*directed_height_a + (base_b + base_c)*directed_height_b)/2;
}

static double randomQuadrilateral(GT_Point points[static 4], gsl_rng *rng){
    const double x = gsl_ran_gaussian(rng, 100), y = gsl_ran_gaussian(rng, 100);
    const double angle = gsl_ran_flat(rng, 0, 2*M_PI);
    const double base_a = gsl_ran_exponential(rng, 10./3), base_b = gsl_ran_exponential(rng, 10./3),
            base_c = gsl_ran_exponential(rng, 10./3);
    const double directed_height_a = gsl_ran_gaussian(rng, 5), directed_height_b = gsl_ran_gaussian(rng, 5);
    return initQuadrilateral(points, (GT_Point){x, y}, angle, base_a, base_b, base_c, directed_height_a,
                             directed_height_b);
}

int main(){
    gsl_rng_env_setup();
    const gsl_rng_type *rng_T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(rng_T);
    size_t passed_tests = 0;
    for(size_t i = 0; i < NUM_TRIALS; ++i){
        GT_Point points[3];
        double a1 = randomTriangle(points, rng);
        double a2 = computeAreaShoestring(3, points);
        double a3 = computeAreaTriangleSine(points);
        double a4 = GT_Point_area(3, points[0], points[1], points[2]);
        double a5 = computeAreaTriangleSine2(points);
        double a6 = computeAreaTriangleProject(points);
        //printf("a1: %f\na5: %f\n", a1, a5);
        if(fabs(a1 - a2) < GT_EPSILON && fabs(a1 - a3) < GT_EPSILON && fabs(a1 - a4) < GT_EPSILON && fabs(a1 - a5) < GT_EPSILON &&
           fabs(a1 - a6) < GT_EPSILON && checkInteriorAngles(3, points)){
            ++passed_tests;
        }
    }
    printf("Passed %zu/%zu triangle tests\n", passed_tests, (size_t)(NUM_TRIALS));
    passed_tests = 0;
    for(size_t i = 0; i < NUM_TRIALS; ++i){
        GT_Point points[4];
        double a = randomQuadrilateral(points, rng);
        double b = computeAreaShoestring(4, points);
        double c = GT_Point_area(4, points[0], points[1], points[2], points[3]);
        if(fabs(a - b) < GT_EPSILON && fabs(a - c) < GT_EPSILON
           && checkInteriorAngles(4, points)){
            ++passed_tests;
        }
    }
    printf("Passed %zu/%zu quadrilateral tests\n", passed_tests, (size_t)(NUM_TRIALS));
    gsl_rng_free(rng);
}

