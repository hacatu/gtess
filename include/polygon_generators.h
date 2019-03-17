#ifndef GTESS_POLYGON_GENERATORS_H
#define GTESS_POLYGON_GENERATORS_H

#define _GNU_SOURCE
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "point.h"
#include "polygon.h"

double initTriangle(GT_Point points[static 3], GT_Point o, double angle, double base_a, double base_b,
                           double directed_height);

double randomTriangle(GT_Point points[static 3], gsl_rng *rng);

double initQuadrilateral(GT_Point points[static 4], GT_Point o, double angle, double base_a, double base_b,
                                double base_c, double directed_height_a, double directed_height_b);
                                
double randomQuadrilateral(GT_Point points[static 4], gsl_rng *rng);

double randomNGonWithDiameter(size_t n, GT_Point points[static n], GT_Point a, GT_Point b, gsl_rng *rng);

#endif

