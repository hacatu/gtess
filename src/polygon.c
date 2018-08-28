#include <stddef.h>
#include <math.h>
#include "point.h"
#include "polygon.h"

int GT_Polygon_checkConvex(size_t n, const GT_Point points[static n], int cw){
	if(n < 3){
		return 1;
	}else for(size_t i = 0; i < n; ++i){
		if((cw ? -1 : 1)*GT_Point_cross_abc(points[(i + 1)%n], points[i], points[(n + i - 1)%n]) < GT_EPSILON){
			return 0;
		}
	}
	return 1;
}

int GT_Polygon_checkOrigin(size_t n, const GT_Point points[static n], int cw){
	if(n < 2){
		return !n || GT_Point_sqmag(points[0]) < GT_EPSILON;
	}else for(size_t i = 0; i < n; ++i){
		if((cw ? -1 : 1)*GT_Point_cross(points[i], points[(i + 1)%n]) < GT_EPSILON){
			return 0;
		}
	}
	return 1;
}

int GT_Polygon_contains(size_t n, const GT_Point points[static n], GT_Point p, int cw){
	if(n < 2){
		return !n || GT_Point_sqdist(points[0], p) < GT_EPSILON;
	}else for(size_t i = 0; i < n; ++i){
		if((cw ? -1 : 1)*GT_Point_cross_abc(points[i], p, points[(i + 1)%n]) < GT_EPSILON){
			return 0;
		}
	}
	return 1;
}

double GT_Polygon_diameter(size_t n, const GT_Point points[static n], int cw){
	if(n < 3){
		switch(n){
		case 0:
			return INFINITY;
		case 1:
			return 0;
		case 2:
			return GT_Point_dist(points[0], points[1]);
		}
	}
	size_t top_i = 2;
	GT_Point bottom_axis = GT_Point_unit(GT_Point_sub(points[1], points[0]));
	GT_Point caliper_axis = GT_Point_ccw(bottom_axis);
	double height = GT_Point_dot(caliper_axis, points[top_i]);
	for(size_t i = 3; i < n; ++i){
		double height_here = GT_Point_dot(caliper_axis, points[i]);
		if(height_here < height){
			break;
		}
		height = height_here;
		top_i = i;
	}
	if(fabs(height) < GT_EPSILON){//the points are colinear
		double min_dot = GT_Point_dot(bottom_axis, points[0]);
		double max_dot = GT_Point_dot(bottom_axis, points[1]);
		for(size_t i = 2; i < n; ++i){
			double dot = GT_Point_dot(bottom_axis, points[i]);
			if(dot < max_dot){
				break;
			}
			max_dot = dot;
		}
		for(size_t i = n - 1; i > 1; --i){
			double dot = GT_Point_dot(bottom_axis, points[i]);
			if(dot > min_dot){
				break;
			}
			min_dot = dot;
		}
		return max_dot - min_dot;
	}
	double sqdiam = GT_Point_sqdist(points[0], points[top_i]);
	for(size_t bottom_i = 1; bottom_i < n; ++bottom_i){
		bottom_axis = GT_Point_unit(GT_Point_sub(points[(bottom_i + 1)%n], points[bottom_i]));
		caliper_axis = GT_Point_ccw(bottom_axis);
		height = GT_Point_dot(caliper_axis, points[top_i]);
		for(size_t i = (top_i + 1)%n; i != bottom_i; i = (i + 1)%n){
			double height_here = GT_Point_dot(caliper_axis, points[i]);
			if(height_here < height){
				break;
			}
			height = height_here;
			top_i = i;
		}
		double sqdiam_here = GT_Point_sqdist(points[bottom_i], points[top_i]);
		if(sqdiam_here > sqdiam){
			sqdiam = sqdiam_here;
		}
	}
	return sqrt(sqdiam);
}

typedef struct IntersectConvexContext IntersectConvexContext;
struct IntersectConvexContext {
	size_t n_a, n_b, n_c;//sizes of polygons
	const GT_Point *points_a;//pointers to counterclockwise point lists
	const GT_Point *points_b;
	GT_Point *points_c;
	size_t a_max_i, b_max_i;//indices of the points in a and b that are farthest to the right
	//size_t ;
	size_t seg_lefts[4];//a_bottom_a, a_top_b, b_bottom_a, b_top_b//indices of the [a, b) ends of the counterclockwise segments currently intersecting the scanline
	size_t seg_rights[4];//a_bottom_b, a_top_a, b_bottom_b, b_top_a//indices of the [a, b) ends of the counterclockwise segments currently intersecting the scanline
	/*enum{
		GT_ICC_bBaA,// a > b in the notes; b_bottom <= b_top <= a_bottom <= a_top
		GT_ICC_baBA,// a_ subset b in the notes; b_bottom <= a_bottom <= b_top <= a_top
		GT_ICC_abBA,// b subset a in the notes; a_bottom <= b_bottom <= b_top <= a_top
		GT_ICC_baAB,// a subset b in the notes; b_bottom <= a_bottom <= a_top <= b_top
		GT_ICC_abAB,// a^ subset b in the notes; a_bottom <= b_bottom <= a_top <= b_top
		GT_ICC_aAbB,// a < b in the notes; a_bottom <= a_top <= b_bottom <= b_top
	} seg_order;*/
	int (*state_fn)(IntersectConvexContext *ctx);
	int found_intersection;
	double scanline_x;//should we add a heap of upcoming points of interest (the ends of the 4 active segments and up to 3 intersections on them)?  This might reduce redundant checks.
};

static inline int intersectConvex_scan_bBaA(IntersectConvexContext *ctx);// a > b in the notes; b_bottom <= b_top <= a_bottom <= a_top
static inline int intersectConvex_scan_baBA(IntersectConvexContext *ctx);// a_ subset b in the notes; b_bottom <= a_bottom <= b_top <= a_top
static inline int intersectConvex_scan_abBA(IntersectConvexContext *ctx);// b subset a in the notes; a_bottom <= b_bottom <= b_top <= a_top
static inline int intersectConvex_scan_baAB(IntersectConvexContext *ctx);// a subset b in the notes; b_bottom <= a_bottom <= a_top <= b_top
static inline int intersectConvex_scan_abAB(IntersectConvexContext *ctx);// a^ subset b in the notes; a_bottom <= b_bottom <= a_top <= b_top
static inline int intersectConvex_scan_aAbB(IntersectConvexContext *ctx);// a < b in the notes; a_bottom <= a_top <= b_bottom <= b_top

static inline int intersectConvex_ensureXOrder(IntersectConvexContext *ctx){
	double a_min_x = ctx->points_a[0].x, a_max_x = a_min_x;
	size_t a_min_i = 0, a_max_i = a_min_i;
	for(size_t i = 1; i < ctx->n_a; ++i){//TODO: this loop appears a lot and should be a function, although we don't always want both the min and max
		double x = ctx->points_a[i].x;//TODO: also it could be taking advantage of the fact that convex polygons are bitonic in any coordinate
		if(x <= a_min_x){//TODO: handling points on the same vertical line is hard because in general colinear points aren't good
			a_min_x = x;
			a_min_i = i;
		}else if(x > a_max_x){
			a_max_x = x;
			a_max_i = i;
		}
	}
	double b_min_x = ctx->points_b[0].x, b_max_x = b_min_x;
	size_t b_min_i = 0, b_max_i = b_min_i;
	for(size_t i = 1; i < ctx->n_b; ++i){
		double x = ctx->points_b[i].x;
		if(x <= b_min_x){
			b_min_x = x;
			b_min_i = i;
		}else if(x > b_max_x){
			b_max_x = x;
			b_max_i = i;
		}
	}
	if(b_min_x < a_min_x){//we need to relabel to make sure b does not start to the left of a
		if(b_max_x < a_min_x){
			return 0;//a and b don't overlap so we don't actually have to swap them
		}
		const GT_Point *tmp_points = ctx->points_a;
		ctx->points_a = ctx->points_b;
		ctx->points_b = tmp_points;
		size_t tmp_n = ctx->n_a;
		ctx->n_a = ctx->n_b;
		ctx->n_b = tmp_n;
		ctx->a_max_i = b_max_i;
		ctx->seg_lefts[0] = ctx->seg_lefts[1] = b_min_i;
		ctx->b_max_i = a_max_i;
		ctx->seg_lefts[2] = ctx->seg_lefts[3] = a_min_i;
	}else if(a_max_x < b_min_x){
		return 0;
	}else{
		ctx->a_max_i = a_max_i;
		ctx->seg_lefts[0] = ctx->seg_lefts[1] = a_min_i;
		ctx->b_max_i = b_max_i;
		ctx->seg_lefts[2] = ctx->seg_lefts[3] = b_min_i;
	}
	ctx->seg_rights[0] = (ctx->seg_lefts[0] + 1)%ctx->n_a;
	ctx->seg_rights[1] = (ctx->seg_lefts[1] + ctx->n_a - 1)%ctx->n_a;
	ctx->seg_rights[2] = (ctx->seg_lefts[2] + 1)%ctx->n_b;
	ctx->seg_rights[3] = (ctx->seg_lefts[3] + ctx->n_b - 1)%ctx->n_b;
	return 1;
}

static inline size_t intersectConvex_flatA(IntersectConvexContext *ctx){
	return 0;//TODO: NYI
}

static inline size_t intersectConvex_adjacentX(IntersectConvexContext *ctx){
	return 0;//TODO: NYI
}

static inline int intersectConvex_scanToBLeft(IntersectConvexContext *ctx){
	double b_left = ctx->points_b[ctx->seg_lefts[3]].x;
	{
		size_t seg_b = ctx->seg_lefts[1];
		size_t seg_a = ctx->seg_rights[1];//this can't go past a_max_i since a and b have been verified to overlap and a has been verified to have nonzero width before calling this function
		double seg_a_x = ctx->points_a[seg_a].x;
		while(seg_a_x < b_left){//guaranteed to terminate correctly, see above note
			seg_b = seg_a;
			seg_a = (seg_b + ctx->n_a - 1)%ctx->n_a;
			seg_a_x = ctx->points_a[seg_a].x;
		}
		ctx->seg_lefts[1] = seg_b;
		ctx->seg_rights[1] = seg_a;
	}
	{
		size_t seg_a = ctx->seg_lefts[0];
		size_t seg_b = ctx->seg_lefts[0];
		double seg_b_x = ctx->points_a[seg_b].x;
		while(seg_b_x < b_left){//in addition to the guarantees for the top loop, we know b_left is not a_right since that is also checked before calling this function so we don't need to handle the seg_b_x == b_left case
			seg_a = seg_b;
			seg_b = (seg_a + 1)%ctx->n_a;
			seg_b_x = ctx->points_a[seg_b].x;
		}
		ctx->seg_lefts[0] = seg_a;
		ctx->seg_rights[0] = seg_b;
	}
	ctx->scanline_x = ctx->points_b[ctx->seg_lefts[3]].x;
	return 1;
}

static inline size_t intersectConvex_scan_findNextSegEnd(IntersectConvexContext *ctx, double *x){
	size_t max_i = 0;
	double max_x = ctx->points_a[ctx->seg_rights[0]].x;
	for(size_t i = 1; i < 4; ++i){
		const GT_Point *points = i < 2 ? ctx->points_a : ctx->points_b;
		double curr_x = points[ctx->seg_rights[i]].x;
		if(curr_x > max_x){
			max_x = curr_x;
			max_i = i;
		}
	}
	*x = max_x;
	return max_i;
}

static inline int intersectConvex_scan_advanceToNextSeg(IntersectConvexContext *ctx, int _i){
	size_t n, max_i;//TODO: this is a cute function but it needs to be able to identify when segments intersect at their endpoints and switch ctx->seg_order accordingly
	const GT_Point *points;
	if(_i < 2){
		n = ctx->n_a;
		points = ctx->points_a;
		max_i = ctx->a_max_i;
	}else{
		n = ctx->n_b;
		points = ctx->points_b;
		max_i = ctx->b_max_i;
	}
	size_t seg_step = _i & 1 ? n - 1 : 1;
	ctx->seg_lefts[_i] = ctx->seg_rights[_i];
	ctx->seg_rights[_i] = (ctx->seg_rights[_i] + seg_step)%n;
	ctx->scanline_x = points[ctx->seg_lefts[_i]].x;
	return ctx->seg_lefts[_i] != max_i;
}

static inline int intersectConvex_scan_checkIntersectionCloser(IntersectConvexContext *ctx, size_t a_seg, size_t b_seg, double *x){
	GT_Point p, a_a, a_b, b_a, b_b;
	if(a_seg == 0){
		a_a = ctx->points_a[ctx->seg_lefts[0]];
		a_b = ctx->points_a[ctx->seg_rights[0]];
	}else{
		a_b = ctx->points_a[ctx->seg_lefts[1]];
		a_a = ctx->points_a[ctx->seg_rights[1]];
	}
	if(b_seg == 2){
		b_a = ctx->points_b[ctx->seg_lefts[2]];
		b_b = ctx->points_a[ctx->seg_rights[2]];
	}else{
		b_b = ctx->points_b[ctx->seg_lefts[3]];
		b_a = ctx->points_b[ctx->seg_rights[3]];
	}
	int status = GT_Point_intersect_segments_ab(&p, a_a, a_b, b_a, b_b);
	if(status){
		if(status == 2){
			a_a = a_seg == 0 ? a_a : a_b;
			b_a = b_seg == 2 ? b_a : b_b;
			p = b_a.x > a_a.x ? b_a : a_a;
			return 2;
		}else if(ctx->scanline_x < p.x && p.x < *x){
			*x = p.x;
			return 1;
		}
	}
	return 0;
}

static inline int intersectConvex_scan_bBaA(IntersectConvexContext *ctx){// a > b in the notes; b_bottom <= b_top <= a_bottom <= a_top
	double x;
	int _i = intersectConvex_scan_findNextSegEnd(ctx, &x);
	//check for intersection between a and B and if its x coordinate is smaller than x then do work with it
	GT_Point p;
	int status = GT_Point_intersect_segments_ab(&p, ctx->points_a[ctx->seg_lefts[0]], ctx->points_a[ctx->seg_rights[0]],
		ctx->points_b[ctx->seg_rights[3]], ctx->points_b[ctx->seg_lefts[3]]);
	if(status){
		if(status == 2){//the segments are colinear so their first real intersection is at ctx->points_a[ctx->seg_lefts[0]].x.  This point is already a vertex, we just have to make sure we are in a central state so it gets added to the output
			ctx->state_fn = intersectConvex_scan_abAB;
			return 1;
		}else if(ctx->scanline_x < p.x && p.x < x){
			ctx->scanline_x = p.x;
			ctx->state_fn = intersectConvex_scan_abAB;
			return 1;
		}else if(p.x == x){
			//TODO: handle intersections at shared vertices or vertices, ie p.x == x
			//for now we will do nothing (since both parts of the intersection might not have been passed) and handle it in advanceToNextSeg.
		}
	}
	return intersectConvex_scan_advanceToNextSeg(ctx, _i);
}

static inline int intersectConvex_scan_baBA(IntersectConvexContext *ctx){// a_ subset b in the notes; b_bottom <= a_bottom <= b_top <= a_top
	double x;
	int _i = intersectConvex_scan_findNextSegEnd(ctx, &x);
	//check for intersections between a and b, A and B, and a and B
	GT_Point p, p1, p2;
	//intersection between a and b going to state b subset a
	int status = GT_Point_intersect_segments_ab(&p, ctx->points_a[ctx->seg_lefts[0]], ctx->points_a[ctx->seg_rights[0]],
		ctx->points_b[ctx->seg_lefts[2]], ctx->points_b[ctx->seg_rights[2]]);
	int (*next_state_fn)(IntersectConvexContext *ctx) = ctx->state_fn;
	if(status){
		if(status == 2){
			p = ctx->points_a[ctx->seg_lefts[0]];
			p1 = ctx->points_b[ctx->seg_lefts[2]];
			if(p1.x > p.x){
				p = p1;
			}
		}
		if(ctx->scanline_x < p.x && p.x < x){
			next_state_fn = intersectConvex_scan_abBA;
		}else if(p.x == x){
			//TODO: handle intersections at vertices
		}
	}
	//intersection of A and B going to state a subset b
	status = GT_Point_intersect_segments_ab(&p1, ctx->points_a[ctx->seg_rights[1]], ctx->points_a[ctx->seg_lefts[1]],
		ctx->points_b[ctx->seg_rights[3]], ctx->points_b[ctx->seg_lefts[3]]);
	if(status){
		if(status == 2){
			p1 = ctx->points_a[ctx->seg_lefts[1]];
			p2 = ctx->points_b[ctx->seg_lefts[3]];
			if(p2.x > p1.x){
				p1 = p2;
			}
		}
		if(ctx->scanline_x < p1.x && p1.x < x){
			if(next_state_fn == ctx->state_fn || p1.x < p.x){
				p = p1;
				next_state_fn = intersectConvex_scan_baAB;
			}
		}else if(p1.x == x){
			//TODO: handle intersections at vertices
		}
	}
	//intersection of a and B going to state bBaA
	status = GT_Point_intersect_segments_ab(&p1, ctx->points_a[ctx->seg_lefts[0]], ctx->points_a[ctx->seg_rights[0]],
		ctx->points_b[ctx->seg_rights[3]], ctx->points_b[ctx->seg_lefts[3]]);
	if(status){
		if(status == 2){
			p1 = ctx->points_a[ctx->seg_lefts[0]];
			p2 = ctx->points_b[ctx->seg_lefts[3]];
			if(p2.x > p1.x){
				p1 = p2;
			}
		}
		if(ctx->scanline_x < p1.x && p1.x < x){
			if(next_state_fn == ctx->state_fn || p1.x < p.x){
				p = p1;
				next_state_fn = intersectConvex_scan_bBaA;
			}
		}else if(p1.x == x){
			//TODO: handle intersections at vertices
		}
	}
	if(next_state_fn == ctx->state_fn){
		return intersectConvex_scan_advanceToNextSeg(ctx, _i);
	}
	ctx->scanline_x = p.x;
	ctx->state_fn = next_state_fn;
	return 1;
}

static inline int intersectConvex_scan_abBA(IntersectConvexContext *ctx){// b subset a in the notes; a_bottom <= b_bottom <= b_top <= a_top
	return 0;//TODO: NYI
}

static inline int intersectConvex_scan_baAB(IntersectConvexContext *ctx){// a subset b in the notes; b_bottom <= a_bottom <= a_top <= b_top
	return 0;//TODO: NYI
}

static inline int intersectConvex_scan_abAB(IntersectConvexContext *ctx){// a^ subset b in the notes; a_bottom <= b_bottom <= a_top <= b_top
	return 0;//TODO: NYI
}

static inline int intersectConvex_scan_aAbB(IntersectConvexContext *ctx){// a < b in the notes; a_bottom <= a_top <= b_bottom <= b_top
	return 0;//TODO: NYI
}

size_t GT_Polygon_intersectConvex(GT_Point *points_c, size_t n_a, const GT_Point points_a[static n_a], size_t n_b, const GT_Point points_b[static n_b]){
	if(!n_a || !n_b){//TODO: Make sure an empyt list of points is always considered as an empty polygon when that makes sense (previously I've used it to mean everything)
		return 0;
	}
	//segment a_top, a_bottom, b_top, b_bottom;
	//we need to find the leftmost points in a and b and swap them so b is not farther left than a
	IntersectConvexContext ctx = {
		.n_a = n_a,
		.points_a = points_a,
		.n_b = n_b,
		.points_b = points_b,
		.points_c = points_c,
	};//nb unspecified members are set to zero (designated initializers)
	if(!intersectConvex_ensureXOrder(&ctx)){
		return 0;
	}else if(ctx.seg_lefts[1] == ctx.a_max_i){//a has no segments from its leftmost point to its rightmost point which means it has zero width//TODO: should these checks have epsilons?  Should all < and <= checks have epsilons?
		return intersectConvex_flatA(&ctx);
	}else if(ctx.points_a[ctx.a_max_i].x == ctx.points_b[ctx.seg_lefts[3]].x){//a and b only overlap at the very end
		return intersectConvex_adjacentX(&ctx);
	}else if(!intersectConvex_scanToBLeft(&ctx)){
		return 0;
	}
	if(ctx.points_a[ctx.seg_lefts[1]].y < ctx.points_b[ctx.seg_lefts[2]].y){
		ctx.state_fn = intersectConvex_scan_aAbB;
		//start case a < b
	}else if(ctx.points_a[ctx.seg_lefts[0]].y > ctx.points_b[ctx.seg_lefts[3]].y){
		ctx.state_fn = intersectConvex_scan_bBaA;
		//start case a > b
	}else{
		ctx.state_fn = intersectConvex_scan_abBA;
		ctx.found_intersection = 1;
		//start case b subset a
	}
	int status = 0;
	do{
		status = ctx.state_fn(&ctx);
	}while(status);
	return ctx.n_c;
}

