#define _GNU_SOURCE
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "point.h"
#include "polygon.h"
#include <stdio.h>

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

double randomConvexNGonWithDiameter(size_t *bi, size_t n, GT_Point points[static n], GT_Point a, GT_Point b, gsl_rng *rng){
	/* algorithm: 2 of the n points are a and b.  The other n - 2 we add must lie within the closed ball
	 * with center (a + b)/2 and radius |b - a|/2.
	 * We can pick n - 2 angles uniformly and then to get points uniformly distributed on the ball we pick the
	 * corresponding radii as the square root of uniformly distributed random variables.
	 * We only use these uniform points on the circle to generate positions along the diameter.
	 * This could also be done by inverting the cdf of a circle using either a series or newton's method (but I
	 * didn't want to do this) or bounding the radii initially by generating cartesian-linear constraints in polar
	 * coordinates.
	 * 
	 * The pdf for a circle with raduis 1 is 2*sqrt(1 - x^2) so it has cdf asin(x) + x*sqrt(1 - x^2) + pi/2 (not normalized).
	 * If we divide by cdf(1) - cdf(-1) = pi we get a normalized cdf F whose taylor series can be found by lagrange inversion.
	 * The even power coefficients will be 0 and the odd terms will be g_n*(x - 1/2)^n/n! where g_1 = pi/2, g_3 = pi^3/8,
	 * g_5 = 13*pi^5/32, g_7 = 493*pi^7/128, g_9 = 37369*pi^9/512, g_11 = 4732249*pi^11/2048
	 * It appears this sequence matches https://oeis.org/A281181
	 * 
	 * I ended up generating n - 2 points on the unit ball and projecting them onto the interval [-1, 1].
	 * These points are sorted by x coordinate.
	 * Then the y coordinates are generated using the line from (-1, 0) through the last top point as an upper bound
	 * and the line from the last top point through (1, 0) as a lower bound for the upper region and similarly
	 * for the lower region.  The number of bottom points is counted as they are created so the top and
	 * bottom points can easily be de-interspersed.  Finally, the actual points are generated by transforming the
	 * points we generated on the unit circle onto the goal circle.  The area of the unit circle polygon is calculated
	 * using the induced trapezoid (and triangle) decomposition and then scaled by r^2.
	 */
	GT_Point c = GT_Point_midpoint(a, b);
	double r = GT_Point_dist(a, b)/2;
	double A = 0, dA;
	for(size_t i = 1; i < n - 1; ++i){
		double x, y;
		do{
			x = gsl_ran_flat(rng, -1., 1.);
			y = gsl_ran_flat(rng, -1., 1.);
		}while(x*x + y*y > 1);
		points[i].x = x;
	}
	qsort(points + 1, n - 2, sizeof(GT_Point), GT_Point_cmp_fn_x);
	points[0] = (GT_Point){-1., 0.};
	points[n - 1] = (GT_Point){1., 0.};
	size_t thead_i = 0, bhead_i = 0, blen = 0;
	double tm = INFINITY, bm = -INFINITY;
	for(size_t i = 1; i < n - 1; ++i){
		double by = fmax(
			bm*(points[i].x - points[bhead_i].x) + points[bhead_i].y,
			-sqrt(1 - points[i].x*points[i].x));
		double Ay = fmin(
			tm*(points[i].x - points[thead_i].x) + points[thead_i].y,
			sqrt(1 - points[i].x*points[i].x));
		double By = points[bhead_i].y*(1 - points[i].x)/(1 - points[bhead_i].x);
		double ay = points[thead_i].y*(1 - points[i].x)/(1 - points[thead_i].x);
		double y = gsl_ran_flat(rng, by, Ay - ay + By);
		if(y > By){
			y += ay - By;
			tm = (y - points[thead_i].y)/(points[i].x - points[thead_i].x);
			dA = (points[thead_i].y + y)*(points[i].x - points[thead_i].x)/2;
			//printf("Area of new top trapezoid: %f\n", dA);
			A += dA;
			thead_i = i;
		}else{
			bm = (y - points[bhead_i].y)/(points[i].x - points[bhead_i].x);
			dA = -(points[bhead_i].y + y)*(points[i].x - points[bhead_i].x)/2;
			//printf("Area of new bottom trapezoid: %f\n", dA);
			A += dA;
			bhead_i = i;
			++blen;
		}
		points[i].y = y;
	}
	dA = points[thead_i].y*(1 - points[thead_i].x)/2;
	//printf("Area of last top trapezoid: %f\n", dA);
	A += dA;
	dA = -points[bhead_i].y*(1 - points[bhead_i].x)/2;
	//printf("Area of last bottom trapezoid: %f\n", dA);
	A += dA;
	GT_Point *bpoints = malloc(blen*sizeof(GT_Point));
	if(!bpoints){
		return -1.;
	}
	//Here I use a clockwise perpendicular vector instead of ccw because the points on the unit ball are arranged clockwise but the image should have them ccw
	GT_Point xd = GT_Point_sub(b, c), yd = GT_Point_cw(xd);
	for(size_t i = 1, j = 1, k = blen - 1; i < n - 1; ++i){
		if(points[i].y > 0){
			points[j++] = GT_Point_add(c, GT_Point_lincomb2(xd, points[i].x, yd, points[i].y));
		}else{
			bpoints[k--] = GT_Point_add(c, GT_Point_lincomb2(xd, points[i].x, yd, points[i].y));
		}
	}
	memcpy(points + n - blen, bpoints, blen*sizeof(GT_Point));
	free(bpoints);
	points[0] = a;
	//b goes at 0 + 1 + tlen.  1 + tlen + 1 + blen = n so tlen = n - blen - 2 so b goes at n - blen - 1
	points[n - 1 - blen] = b;
	if(bi){
		*bi = n - 1 - blen;
	}
	return A*r*r;
}

double randomConvexNGon(size_t *bi, size_t n, GT_Point points[static n], gsl_rng *rng){
	const double x = gsl_ran_gaussian(rng, 100), y = gsl_ran_gaussian(rng, 100);
	const double dx = gsl_ran_gaussian(rng, 10), dy = gsl_ran_gaussian(rng, 10);
	const GT_Point a = {x, y}, d = {dx, dy}, b = GT_Point_add(a, d);
	return randomConvexNGonWithDiameter(bi, n, points, a, b, rng);
}

GT_Point randomPointIn(size_t n, size_t bi, const GT_Point points[static n], double A, gsl_rng *rng){
	/* we pretty much want to generate a list of trapezoids/triangles and areas similar to how randomConvexNGonWithDiameter finds the area.
	 * Since this function takes the area as a parameter, we can generate a random number between 0 and A and pick the first trapezoid so the
	 * sum of the previous trapezoid's areas is less than the random number but the sum plus the current area.
	 * Then we can pick a random point in that trapezoid pretty easily.
	 */
	double w = gsl_ran_flat(rng, 0, A), dA;
	//for now the 2 or 4 triangles we get (the first and last trapezoids for the top and bottom) will not be special cased, there is a minimal benefit to doing so
	GT_Point a = points[0], d = GT_Point_unit(GT_Point_sub(points[bi], a)), dt = GT_Point_cw(d);
	for(size_t i = 0; i + 1 < n; ++i){
		GT_Point pi = GT_Point_sub(points[i], a);
		double pil = GT_Point_dot(pi, d);//p_i parallel part
		double pit = GT_Point_dot(pi, dt);//p_i perpendicular part
		GT_Point pi1 = GT_Point_sub(points[i + 1], a);
		double pi1l = GT_Point_dot(pi1, d);
		double pi1t = GT_Point_dot(pi1, dt);
		double b = pi1l - pil;
		dA = (pit + pi1t)*b/2;
		if(w < dA){
			w /= dA;//since a uniform(0, A) random variable is equivalent to picking a trapezoid weighted by area dA and then a uniform(0, dA) rv,
			//we normalize the uniform(0, dA) part to a uniform(0, 1) part.  Then if y is the dt coordinate and x the d coordinant of the point,
			//we can pick y uniformly from the slice of the trapezoid at d=x once we have y, so we have to compute the marginal cdf for x and invert
			//it so we can turn a probability index (ie uniform(0, 1) rv such as w) into an evenly distributed x value.
			//if a is the trailing side length and c is the leading side length,
			//P(d <= x) = (area of partial trapezoid up to x)/(area of entire trapezoid) = (a + ((c - a)*x/b + a))*x/2 / (dA = (a + c)*b/2)
			//= (2*b*a*x + (c - a)*x^2 / (a + c)*b^2, and we invert to find x in terms of P,
			//0 = (c - a)*x^2 + 2*b*a*x - (a + c)*b^2*P, x = (-2*b*a +/- sqrt(4*b^2*a^2 + 4*b^2*(c - a)*(a + c)*P))/2*(c - a)
			//= b*(-a +/- sqrt(a^2 + (c - a)*(a + c)*P))/(c - a).  So now we just have to pick + or - (also the c = a case should be accounted for)
			//dt is a clockwise transpose of d so this is a left handed coordinate basis (I'm not sure why I picked that instead of a ccw rotation,
			//it was probably because in graphics the y axis is left handed) and we start going from points[0] counterclockwise around the polygon through points[bi]
			//and all the way to points[n-1].
			//This means for i + 1 <= b we have trapezoids under the diameter so the d coordinate of the sides is increasing and they are on the positive dt side.
			//For i >= b we have trapezoids over the diameter so the d coordinate of the sides is decreasing and they are on the negative dt side.
			//This means b = pi1l - pil, the change in the d coordiate from point i to point i + 1, has the same sign as pit and pi1t.
			//We need x to have the same sign as b = pi1l - pil so that when we move x in the d direction we move from point i towards point i + 1.
			//When a and c are positive (the before bi case), assume c > a so c^2 > a^2 so x is positive if we pick +.
			//When a and c are negative (the after bi case), assume c < a so c^2 > a^2 so the numerator (besides b) has the sign we choose so x is negative if we pick -.
			//I checked those and the flipped a vs c cases by graphing the x vs P curve and it checks out.
			double x = b*(copysign(sqrt(pi1t*pi1t + w*(pit - pi1t)*(pit + pi1t)), b) - pi1t)/(pit - pi1t);
			//our randomly selected point will have d coordinate x starting from pi projected onto the diameter
			//now the "top" of the trapezoid at x is pit + (pi1t - pit)*x/b
			double y = (pit + (pi1t - pit)*x/b)*gsl_ran_flat(rng, 0, 1);
			return GT_Point_add(a, GT_Point_lincomb2(d, x + pil, dt, y));
		}
		w -= dA;
	}
	//if we have reached the last trapezoid, it is actually one of the triangles (the first trapezoid is also a triangle but I didn't special case it).
	//We compuete more simply as an offset from a
	GT_Point pi = GT_Point_sub(points[n - 1], a);
	double pil = GT_Point_dot(pi, d);//p_i parallel part
	double pit = GT_Point_dot(pi, dt);//p_i perpendicular part
	dA = -pit*pil/2;
	w /= dA;
	double x = pil*sqrt(w);
	//our randomly selected point will have d coordinate x starting from a
	double y = x*pit*gsl_ran_flat(rng, 0, 1)/pil;
	return GT_Point_add(a, GT_Point_lincomb2(d, x, dt, y));
}

GT_Point randomPointGaussian(GT_Point o, double r, gsl_rng *rng){
	return (GT_Point){.x= o.x + gsl_ran_gaussian(rng, r), .y= o.y + gsl_ran_gaussian(rng, r)};
}

