#define _GNU_SOURCE
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "point.h"

#define NUM_TRIALS 100

int main(){
    gsl_rng_env_setup();
    const gsl_rng_type *rng_T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(rng_T);
    size_t passed_tests = 0;
    for(size_t i = 0; i < NUM_TRIALS; ++i){
        double t = gsl_ran_flat(rng, 0, 2*M_PI);
        double cos_t = cos(t), sin_t = sin(t);
        GT_Point p = GT_Point_lincomb2(GT_Point_e2, sin_t, GT_Point_e1, cos_t);
        if(fabs(GT_Point_mag(p) - 1) >= GT_EPSILON){
            continue;
        }else if(GT_Point_sqdist(GT_Point_lincomb(2, GT_Point_e1, cos_t, GT_Point_e2, sin_t), p) >= GT_EPSILON){
            continue;
        }else if(GT_Point_sqdist(GT_Point_scale(GT_Point_e1, cos_t), GT_Point_project_unit(p, GT_Point_e1)) >= GT_EPSILON){
            continue;
        }
        ++passed_tests;
    }
    printf("Passed %zu/%zu tests\n", passed_tests, (size_t)(NUM_TRIALS));
    gsl_rng_free(rng);
}

