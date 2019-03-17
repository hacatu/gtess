#define _GNU_SOURCE
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "point.h"
#include "polygon.h"

double initTriangle(GT_Point points[static 3], GT_Point o, double angle, double base_a, double base_b,
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

double randomTriangle(GT_Point points[static 3], gsl_rng *rng){
    const double x = gsl_ran_gaussian(rng, 100), y = gsl_ran_gaussian(rng, 100);
    const double angle = gsl_ran_flat(rng, 0, 2*M_PI);
    const double base_a = gsl_ran_exponential(rng, 5), base_b = gsl_ran_exponential(rng, 5);
    const double directed_height = gsl_ran_gaussian(rng, 10);
    return initTriangle(points, (GT_Point){x, y}, angle, base_a, base_b, directed_height);
}

double initQuadrilateral(GT_Point points[static 4], GT_Point o, double angle, double base_a, double base_b,
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

double randomQuadrilateral(GT_Point points[static 4], gsl_rng *rng){
    const double x = gsl_ran_gaussian(rng, 100), y = gsl_ran_gaussian(rng, 100);
    const double angle = gsl_ran_flat(rng, 0, 2*M_PI);
    const double base_a = gsl_ran_exponential(rng, 10./3), base_b = gsl_ran_exponential(rng, 10./3),
            base_c = gsl_ran_exponential(rng, 10./3);
    const double directed_height_a = gsl_ran_gaussian(rng, 5), directed_height_b = gsl_ran_gaussian(rng, 5);
    return initQuadrilateral(points, (GT_Point){x, y}, angle, base_a, base_b, base_c, directed_height_a,
                             directed_height_b);
}

double randomNGonWithDiameter(size_t n, GT_Point points[static n], GT_Point a, GT_Point b, gsl_rng *rng){
	/* algorithm: 2 of the n points are a and b.  The other n - 2 we add must lie within the closed ball
	 * with center (a + b)/2 and radius |b - a|/2.
	 * We can pick n - 2 angles uniformly, rejecting any which are too close.
	 */
	//TODO: deal with n < 2 and a == b
	double angle;
	for(size_t i = 2; i < n; ++i){
		do{
			angle = gsl_ran_flat(rng, 0, 2*M_PI);
		}while(fmod(angle, M_PI) >= GT_EPSILON);
		points[i].x = angle;
	}
	for(int dups = 1; dups;){
		dups = 0;
		qsort(points + 2, n - 2, sizeof(GT_Point), GT_Point_cmp_fn_x);
		for(size_t i = 3; i < n; ++i){
			if(points[i].x - points[i - 1].x < GT_EPSILON){
				dups = 1;
				do{
					angle = gsl_ran_flat(rng, 0, 2*M_PI);
				}while(fmod(angle, M_PI) >= GT_EPSILON);
				points[i].x = angle;
			}
		}
	}
	points[0] = GT_Point_e2;
	points[1] = GT_Point_neg1_polar;
	qsort(points, n, sizeof(GT_Point), GT_Point_cmp_fn_x);
	//Now we have n unique angles, 2 of which are 0 and pi corresponding to a and b respectively, in sorted order.
	//For the first and last, we can pick any radius.  Then, we find the index of pi and work backwards.
	//In general, if we know 
	return -1.;
}

