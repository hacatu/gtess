#define _GNU_SOURCE
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "point.h"
#include "polygon.h"
#include "polygon_generators.h"

#define NUM_TRIALS 100

int main(){
    gsl_rng_env_setup();
    const gsl_rng_type *rng_T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(rng_T);
    GT_Point points[] = {{0, 0}, {3, -1}, {4, 3}, {2, 3}, {0.5, 1}};//this is ccw and the diameter is (0, 0), (4, 3)
    double A = GT_Polygon_area(5, points);
    GT_Point p = randomPointIn(5, 2, points, A, rng);
    printf("(%f, %f)\n", p.x, p.y);
    /*
    size_t passed_tests = 0;
    for(size_t i = 0; i < NUM_TRIALS; ++i){
        GT_Point points[3];
        double a1 = randomTriangle(points, rng);
        double a2 = GT_Polygon_area(3, points);
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
        double b = GT_Polygon_area(4, points);
        double c = GT_Point_area(4, points[0], points[1], points[2], points[3]);
        if(fabs(a - b) < GT_EPSILON && fabs(a - c) < GT_EPSILON
           && checkInteriorAngles(4, points)){
            ++passed_tests;
        }
    }
    printf("Passed %zu/%zu quadrilateral tests\n", passed_tests, (size_t)(NUM_TRIALS));
    */
    gsl_rng_free(rng);
}

