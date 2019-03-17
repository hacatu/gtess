#define _GNU_SOURCE
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "point.h"
#include "polygon.h"
#include "polygon_generators.h"

int main(){
    gsl_rng_env_setup();
    const gsl_rng_type *rng_T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(rng_T);
    
    gsl_rng_free(rng);
}

