#include <stddef.h>
#include <string.h>
#include <math.h>
#include "point.h"
#include "polygon.h"

typedef struct{
	GT_Point *buf;
	size_t len, a, b, cap;
} GT_Point_Queue;

typedef struct GT_SLE_State GT_SLE_State;
typedef struct GT_SLE_Event GT_SLE_Event;

typedef struct{
	const GT_Point *points;
	size_t n;
	size_t seg_ends[2][2];//[bottom,top][left,right]
	size_t max_i;
} GT_SLE_Polygon;

typedef struct{
	GT_SLE_Event *buf;
	size_t len, cap;
} GT_EventHeap;

struct GT_SLE_State{
	GT_SLE_Polygon polygons[2];
	GT_Point_Queue points;
	size_t boundary_intersections;
	GT_EventHeap events;
	int (*state_fn)(GT_SLE_State *state);
};

struct GT_SLE_Event{
	GT_Point p;
	enum{
		GT_SLE_a_END = 1,
		GT_SLE_A_END = 2,
		GT_SLE_b_END = 4,
		GT_SLE_B_END = 8,
		GT_SLE_ab_CROSS = 16,
		GT_SLE_aB_CROSS = 32,
		GT_SLE_Ab_CROSS = 64,
		GT_SLE_AB_CROSS = 128,
	} conditions;
};

static const int GT_SLE_ANY_END = GT_SLE_a_END | GT_SLE_A_END | GT_SLE_b_END | GT_SLE_B_END;
static const int GT_SLE_ANY_CROSS = GT_SLE_ab_CROSS | GT_SLE_aB_CROSS | GT_SLE_Ab_CROSS | GT_SLE_AB_CROSS;
static const int GT_SLE_CROSS_table[4][4] = {
	[0][2] = GT_SLE_ab_CROSS,
	[2][0] = GT_SLE_ab_CROSS,
	[0][3] = GT_SLE_aB_CROSS,
	[3][0] = GT_SLE_aB_CROSS,
	[1][2] = GT_SLE_Ab_CROSS,
	[2][1] = GT_SLE_Ab_CROSS,
	[1][3] = GT_SLE_AB_CROSS,
	[3][1] = GT_SLE_AB_CROSS,
};

static inline int GT_Point_Queue_pusha(GT_Point_Queue *self, GT_Point p){
	if(self->len == self->cap){
		return 0;
	}
	self->a = (self->a + self->cap - 1)%self->cap;
	self->buf[self->a] = p;
	++self->len;
	return 1;
}

static inline int GT_Point_Queue_pushb(GT_Point_Queue *self, GT_Point p){
	if(self->len == self->cap){
		return 0;
	}
	self->buf[self->b] = p;
	self->b = (self->b + 1)%self->cap;
	++self->len;
	return 1;
}

static inline void GT_Point_Queue_canonicalize(GT_Point_Queue *self){
	if(!self->len || !self->a){
		return;
	}else if(self->a < self->b){
		memmove(self->buf, self->buf + self->a, self->len*sizeof(GT_Point));
	}else if(self->cap - self->a + self->b < self->a + 1){
		memmove(self->buf + self->cap - self->a, self->buf, self->b*sizeof(GT_Point));
		memcpy(self->buf, self->buf + self->a, (self->cap - self->a)*sizeof(GT_Point));
	}else{
		size_t orbits = GT_gcd(self->cap, self->cap - self->a);
		for(size_t i = 0, j = i; i < orbits; j = ++i){
			GT_Point p = self->buf[i];
			while(1){
				size_t j_next = (j + self->cap - self->a)%self->cap;
				self->buf[j] = self->buf[j_next];
				if(j_next == i){
					break;
				}
				j = j_next;
			}
			self->buf[j] = p;
		}
	}
	self->a = 0;
	self->b = self->len;
}

static inline void GT_EventHeap_swap(GT_EventHeap *self, size_t i, size_t j){
	GT_SLE_Event t = self->buf[i];
	self->buf[i] = self->buf[j];
	self->buf[j] = t;
}

static inline void GT_EventHeap_siftup(GT_EventHeap *self, size_t i){
	for(size_t j = (i - 1) >> 1; i; i = j, j = (i - 1) >> 1){
		if(GT_Point_cmp_xy(self->buf[i].p, self->buf[j].p) >= 0){
			break;
		}
		GT_EventHeap_swap(self, i, j);
	}
}

static inline void GT_EventHeap_siftdown(GT_EventHeap *self, size_t i){
	size_t l = 2*i + 1, r = l + 1;
	for(size_t j; r < self->len; i = j, l = 2*i + 1, r = l + 1){
		j = GT_Point_cmp_xy(self->buf[l].p, self->buf[r].p) < 0 ? l : r;
		if(GT_Point_cmp_xy(self->buf[j].p, self->buf[i].p) >= 0){
			break;
		}
		GT_EventHeap_swap(self, i, j);
	}
	if(l < self->len && GT_Point_cmp_xy(self->buf[l].p, self->buf[i].p) < 0){
		GT_EventHeap_swap(self, i, l);
	}
}

static inline int GT_EventHeap_push(GT_EventHeap *self, GT_SLE_Event ev){
	if(self->len == self->cap){
		return 0;
	}
	self->buf[self->len] = ev;
	GT_EventHeap_siftup(self, self->len++);
	return 1;
}

static inline int GT_EventHeap_pop(GT_EventHeap *self, GT_SLE_Event *out){
	if(!self->len){
		return 0;
	}
	if(out){
		*out = self->buf[0];
	}
	if(--self->len){
		self->buf[0] = self->buf[self->len];
		GT_EventHeap_siftdown(self, 0);
	}
	return 1;
}

static inline size_t GT_EventHeap_linsearch(GT_EventHeap *self, GT_Point p){
	for(size_t i = 0; i < self->len; ++i){
		GT_Point q = self->buf[i].p;
		if(q.x == p.x && q.y == p.y){
			return i;
		}
	}
	return self->len;
}

static inline int GT_SLE_iterate_seg(GT_SLE_State *state, int which){
	int is_b = !(which & (GT_SLE_a_END | GT_SLE_A_END));
	GT_SLE_Polygon *polygon = state->polygons + is_b;
	GT_SLE_Polygon *other = state->polygons + !is_b;
	int is_top = !(which & (GT_SLE_a_END | GT_SLE_b_END));
	size_t step = is_top ? polygon->n - 1 : 1;
	size_t *seg_ends = polygon->seg_ends[is_top];
	seg_ends[0] = seg_ends[1];
	if(seg_ends[0] == polygon->max_i){
		return 0;
	}
	seg_ends[1] = (seg_ends[1] + step)%polygon->n;
	GT_Point p = polygon->points[seg_ends[1]];
	size_t i = GT_EventHeap_linsearch(&state->events, p);
	if(i == state->events.len){
		GT_EventHeap_push(&state->events, (GT_SLE_Event){.p=p, .conditions=which});
	}else{
		state->events.buf[i].conditions |= which;
	}
	for(size_t j = 0; j < 2; ++j){
		size_t *other_ends = other->seg_ends[j];
		if(GT_Point_intersect_segments_open_ab(&p, polygon->points[seg_ends[0]], polygon->points[seg_ends[1]],
		                                       other->points[other_ends[0]], other->points[other_ends[1]]) == 1){
			GT_EventHeap_push(&state->events, (GT_SLE_Event){.p=p, .conditions= GT_SLE_CROSS_table[2*is_b + is_top][2*!is_b + j]});
		}//we don't need to check for overlap anymore because overlap of an intersection of open segments implies the end of one of the polygons since they are convex and nondegenerate
	}
	return 1;
}

static inline int GT_SLE_init_events(GT_SLE_State *state){
	GT_Point p;
	for(int which = 1; which & GT_SLE_ANY_END; which <<= 1){
		int is_b = !(which & (GT_SLE_a_END | GT_SLE_A_END));
		int is_top = !(which & (GT_SLE_a_END | GT_SLE_b_END));
		GT_SLE_Polygon *polygon = state->polygons + is_b;
		p = polygon->points[polygon->seg_ends[is_top][1]];
		size_t i = GT_EventHeap_linsearch(&state->events, p);
		if(i == state->events.len){
			GT_EventHeap_push(&state->events, (GT_SLE_Event){.p=p, .conditions=which});
		}else{
			state->events.buf[i].conditions |= which;
		}
	}
	for(int a_top = 0; a_top < 2; ++a_top){
		size_t *a_ends = state->polygons[0].seg_ends[a_top];
		for(int b_top = 0; b_top < 2; ++b_top){
			size_t *b_ends = state->polygons[1].seg_ends[b_top];
			if(GT_Point_intersect_segments_open_ab(&p, state->polygons[0].points[a_ends[0]], state->polygons[0].points[a_ends[1]],
			                                       state->polygons[1].points[b_ends[0]], state->polygons[1].points[b_ends[1]]) == 1){
				GT_EventHeap_push(&state->events, (GT_SLE_Event){.p=p, .conditions= GT_SLE_CROSS_table[a_top][2 + b_top]});
			}
		}
	}
	return 1;
}

static inline int GT_SLE_scan_stateless(GT_SLE_State *state){
	GT_SLE_Event e;
	if(!GT_EventHeap_pop(&state->events, &e)){
		return 0;
	}
	double seg_ys[2][2];
	for(size_t i = 0; i < 2; ++i){
		const GT_SLE_Polygon *polygon = &state->polygons[i];
		for(size_t j = 0; j < 2; ++j){
			GT_Point a = polygon->points[polygon->seg_ends[j][j]];
			GT_Point b = polygon->points[polygon->seg_ends[j][!j]];
			double m = GT_Point_slope(a, b);
			if(isnan(m)){
				seg_ys[i][j] = ((b.y > a.y) == j) ? b.y : a.y;//pick the higher point if we are on top (j == 1) and the lower point if we are on bottom
			}else{
				seg_ys[i][j] = a.y + m*(e.p.x - a.x);
			}
		}
	}
	double max_overlap = seg_ys[0][1] > seg_ys[1][1] ? seg_ys[1][1] : seg_ys[0][1];
	double min_overlap = seg_ys[0][0] < seg_ys[1][0] ? seg_ys[1][0] : seg_ys[0][0];
	if(min_overlap - max_overlap <= GT_EPSILON){//min_overlap <= max_overlap but accounting for error in a way where GT_EPSILON doesn't get absorbed
		if(fabs(e.p.y - min_overlap) < GT_EPSILON){
			GT_Point_Queue_pushb(&state->points, e.p);
		}else if(fabs(e.p.y - max_overlap) < GT_EPSILON){
			GT_Point_Queue_pusha(&state->points, e.p);
		}
	}else if(state->points.len){
		return 0;
	}
	if(e.conditions & GT_SLE_ANY_CROSS){
		++state->boundary_intersections;
	}
	int status = 1;
	for(int which = 1; which & GT_SLE_ANY_END; which <<= 1){
		if(which & e.conditions){
			status &= GT_SLE_iterate_seg(state, which);
		}
	}
	return status;
}

static inline int GT_SLE_prepare_polygon(GT_SLE_Polygon *polygon){
	size_t min_i = 0, max_i = 0;
	GT_Point min_xy = polygon->points[0], max_xy = min_xy;
	for(size_t i = 1; i < polygon->n; ++i){
		GT_Point p = polygon->points[i];
		if(GT_Point_cmp_xy(p, min_xy) < 0){
			min_xy = p;
			min_i = i;
		}else if(GT_Point_cmp_xy(p, max_xy) > 0){
			max_xy = p;
			max_i = i;
		}
	}
	polygon->seg_ends[0][0] = polygon->seg_ends[1][0] = min_i;
	polygon->seg_ends[0][1] = (min_i + 1)%polygon->n;
	polygon->seg_ends[1][1] = (min_i + polygon->n - 1)%polygon->n;
	polygon->max_i = max_i;
	return GT_Polygon_dimension(polygon->n, polygon->points);
}

static inline int GT_SLE_scan_to_start(GT_SLE_State *state, int top){
	GT_SLE_Polygon *polygon = state->polygons;
	size_t step = top ? polygon->n - 1 : 1;
	double b_left = state->polygons[1].points[state->polygons[1].seg_ends[1][0]].x;
	while(polygon->points[polygon->seg_ends[top][1]].x < b_left){
		if(polygon->seg_ends[top][1] == polygon->max_i){
			return 0;
		}
		polygon->seg_ends[top][0] = polygon->seg_ends[top][1];
		polygon->seg_ends[top][1] = (polygon->seg_ends[top][1] + step)%polygon->n;
	}
	return 1;
}

static inline void swapPolygons(GT_SLE_State *state){
	GT_SLE_Polygon t = state->polygons[0];
	state->polygons[0] = state->polygons[1];
	state->polygons[1] = t;
}

static inline int intersectPointPoint(GT_SLE_State *state){
	if(!state->polygons[0].n || !state->polygons[1].n){
		return 0;
	}else if(!GT_Point_cmp_xy(state->polygons[0].points[0], state->polygons[1].points[0])){
		GT_Point_Queue_pushb(&state->points, state->polygons[0].points[0]);
		return 1;
	}
	return 0;
}

static inline int intersectPointLine(GT_SLE_State *state){
	if(!state->polygons[0].n){
		return 0;
	}else if(GT_Point_cross_abc(state->polygons[0].points[0], state->polygons[1].points[state->polygons[1].seg_ends[1][0]], state->polygons[1].points[state->polygons[1].max_i])){
		return 0;
	}else if(GT_Point_cmp_xy(state->polygons[0].points[0], state->polygons[1].points[state->polygons[1].seg_ends[1][0]]) < 0){
		return 0;
	}else if(GT_Point_cmp_xy(state->polygons[0].points[0], state->polygons[1].points[state->polygons[1].max_i]) > 0){
		return 0;
	}
	GT_Point_Queue_pushb(&state->points, state->polygons[0].points[0]);
	return 1;
}

static inline int intersectLinePoint(GT_SLE_State *state){
	swapPolygons(state);
	return intersectPointLine(state);
}

static inline int intersectLineLine(GT_SLE_State *state){
	GT_Point p;
	switch(GT_Point_intersect_segments_closed_ab(&p,
		state->polygons[0].points[state->polygons[0].seg_ends[1][0]],
		state->polygons[0].points[state->polygons[0].max_i],
		state->polygons[1].points[state->polygons[1].seg_ends[1][0]],
		state->polygons[1].points[state->polygons[1].max_i])){
		case 0: return 0;
		case 1:
			GT_Point_Queue_pushb(&state->points, p);
			state->boundary_intersections = 1;
			return 1;
		case 2:
		{
			GT_Point ext_xy[2] = {state->polygons[0].points[state->polygons[0].seg_ends[1][0]],
				state->polygons[0].points[state->polygons[0].max_i]};
			p = state->polygons[1].points[state->polygons[1].seg_ends[1][0]];
			if(GT_Point_cmp_xy(ext_xy[0], p) < 0){
				ext_xy[0] = p;
			}
			p = state->polygons[1].points[state->polygons[1].max_i];
			if(GT_Point_cmp_xy(ext_xy[1], p) > 0){
				ext_xy[1] = p;
			}
			int ord = GT_Point_cmp_xy(ext_xy[0], ext_xy[1]);
			for(size_t i = 0; i < 1 - ord; ++i){//append 0 points if there is no overlap, 1 point if the overlap is 1 point, and 2 points otherwise (ord can be -1, 0, or +1)
				GT_Point_Queue_pushb(&state->points, ext_xy[i]);
			}
			state->boundary_intersections = 2;
			return 1;
		}
		default: __builtin_unreachable();
	}
}

static inline int intersectLineConvex(GT_SLE_State *state){
	GT_Point crosses[2];
	size_t crosses_len = 0;
	GT_Point line_ends[2] = {state->polygons[0].points[state->polygons[0].seg_ends[1][0]],
		state->polygons[0].points[state->polygons[0].max_i]};
	const GT_Point *points = state->polygons[1].points;
	size_t n = state->polygons[1].n, off_i = 0;
	double last_cross;
	while(1){
		if((last_cross = GT_Point_cross_abc(points[off_i], line_ends[0], line_ends[1]))){
			break;
		}
		++off_i;
	};
	for(size_t i = off_i + 1, j; i != off_i; i = (i + 1)%n){
		double cross = GT_Point_cross_abc(points[i], line_ends[0], line_ends[1]);
		if(cross == 0){
			if(crosses_len){
				crosses[1] = points[i];
				crosses_len = 2;
				break;
			}else{
				crosses[0] = points[i];
				crosses_len = 1;
			}
			for(j = (i + 1)%n; j != off_i; j = (j + 1)%n){
				cross = GT_Point_cross_abc(points[j], line_ends[0], line_ends[1]);
				if(cross == 0){
					crosses[1] = points[j];
					crosses_len = 2;
				}else if(j != (i + 1)%n){
					if(GT_Point_cmp_xy(crosses[0], crosses[1]) > 0){
						GT_Point tmp = crosses[0];
						crosses[0] = crosses[1];
						crosses[1] = tmp;
					}
					if(GT_Point_cmp_xy(crosses[0], line_ends[0]) < 0){
						crosses[0] = line_ends[0];
					}
					if(GT_Point_cmp_xy(crosses[1], line_ends[1]) > 0){
						crosses[1] = line_ends[1];
					}
					goto FOUND_TWO_INTERSECTIONS;
				}else{
					i = j;
					break;
				}
			}
		}else if(signbit(cross) != signbit(last_cross)){
			int status = GT_Point_intersect_segments_closed_ab(crosses + crosses_len, line_ends[0], line_ends[1], points[(i + n - 1)%n], points[i]);
			if(status){
				if(++crosses_len == 2){
					break;
				}
			}
		}
		last_cross = cross;
	}
	if(crosses_len == 2 && (GT_Point_cmp_xy(crosses[1], line_ends[0]) < 0 || GT_Point_cmp_xy(crosses[1], line_ends[1]) > 0)){
		--crosses_len;
	}
	if(crosses_len && (GT_Point_cmp_xy(crosses[0], line_ends[0]) < 0 || GT_Point_cmp_xy(crosses[0], line_ends[1]) > 0)){
		crosses[0] = crosses[1];
		--crosses_len;
	}
	state->boundary_intersections = crosses_len;
	switch(crosses_len){
		case 0:
			if(GT_Polygon_contains(n, points, line_ends[0], 0)){
				for(size_t i = 0; i < 2; ++i){
					GT_Point_Queue_pushb(&state->points, line_ends[i]);
				}
				return 1;
			}
			return 0;
		case 1:
			GT_Point_Queue_pushb(&state->points, crosses[0]);
			for(size_t i = 0; i < 2; ++i){
				if(GT_Polygon_contains(n, points, line_ends[i], 0)){
					GT_Point_Queue_pushb(&state->points, line_ends[i]);
					return 1;
				}
			}
			return 1;
		case 2:
			FOUND_TWO_INTERSECTIONS:;
			for(size_t i = 0; i < 2; ++i){
				GT_Point_Queue_pushb(&state->points, crosses[i]);
			}
			return 1;
		default: __builtin_unreachable();
	}
}

static inline int intersectConvexLine(GT_SLE_State *state){
	swapPolygons(state);
	return intersectLineConvex(state);
}

static inline int intersectConvexConvex(GT_SLE_State *state){
	if(state->polygons[1].points[state->polygons[1].seg_ends[1][0]].x < state->polygons[0].points[state->polygons[0].seg_ends[1][0]].x){
		swapPolygons(state);
	}
	if(!GT_SLE_scan_to_start(state, 0)){
		return 0;
	}else if(!GT_SLE_scan_to_start(state, 1)){
		return 0;
	}
	GT_SLE_init_events(state);
	while(state->state_fn(state));
	return !!state->points.len;
}

static inline int intersectPointConvex(GT_SLE_State *state){
	if(!state->polygons[0].n){
		return 0;
	}else if(!GT_Polygon_contains(state->polygons[1].n, state->polygons[1].points, state->polygons[0].points[0], 0)){
		return 0;
	}
	GT_Point_Queue_pushb(&state->points, state->polygons[0].points[0]);
	return 1;
}

static inline int intersectConvexPoint(GT_SLE_State *state){
	swapPolygons(state);
	return intersectPointConvex(state);
}

static int (*const intersectors[3][3])(GT_SLE_State *state) = {
	[0][0]=intersectPointPoint,
	[0][1]=intersectPointLine,
	[1][0]=intersectLinePoint,
	[0][2]=intersectPointConvex,
	[2][0]=intersectConvexPoint,
	[1][1]=intersectLineLine,
	[1][2]=intersectLineConvex,
	[2][1]=intersectConvexLine,
	[2][2]=intersectConvexConvex,
};

size_t GT_Polygon_intersectConvex(GT_Point *points_c, size_t *boundary_intersections,
                                  size_t n_a, const GT_Point points_a[static n_a], size_t n_b, const GT_Point points_b[static n_b]){
	GT_SLE_Event _event_buf[8];
	GT_SLE_State state = {
		.polygons[0]={.points=points_a, .n=n_a},
		.polygons[1]={.points=points_b, .n=n_b},
		.points={.buf=points_c, .cap= n_a + n_b},
		.boundary_intersections=0,
		.events={.buf=_event_buf, .cap=8},
		.state_fn=GT_SLE_scan_stateless,
	};
	int a_dim = GT_SLE_prepare_polygon(state.polygons);
	int b_dim = GT_SLE_prepare_polygon(state.polygons + 1);
	intersectors[a_dim][b_dim](&state);
	GT_Point_Queue_canonicalize(&state.points);
	if(boundary_intersections){
		*boundary_intersections = state.boundary_intersections;
	}
	return state.points.len;
}

int GT_Halfplane_intersectHalfplane(GT_Point *out, GT_Halfplane a, GT_Halfplane b){
	GT_Point a_p = a.point, b_p = b.point, a_o = GT_Point_ccw(a.in_normal), b_o = GT_Point_ccw(b.in_normal);
	double det = GT_Point_cross(a_o, b_o);
	if(fabs(det) < GT_EPSILON){
		if(fabs(GT_Point_cross(GT_Point_sub(b_p, a_p), a_o)) < GT_EPSILON){//collinear
			*out = GT_Point_midpoint(a_p, b_p);
			return 2;
		}
		return 0;
	}
	double t = GT_Point_cross(GT_Point_sub(b_p, a_p), a_o)/det;
	*out = GT_Point_add(a_p, GT_Point_scale(a_o, t));
	return 1;
}

int GT_Halfplane_intersectSegment(GT_Point *out, GT_Halfplane halfplane, GT_Point a, GT_Point b){
	GT_Point a_p = halfplane.point, b_p = a, a_o = GT_Point_ccw(halfplane.in_normal), b_o = GT_Point_sub(b, a);
	double det = GT_Point_cross(a_o, b_o);
	if(fabs(det) < GT_EPSILON){
		if(fabs(GT_Point_cross(GT_Point_sub(b_p, a_p), a_o)) < GT_EPSILON){//collinear
			*out = GT_Point_midpoint(a, b);
			return 2;
		}
		return 0;
	}
	double t = GT_Point_cross(GT_Point_sub(b_p, a_p), a_o)/det;
	double u = GT_Point_cross(GT_Point_sub(b_p, a_p), b_o)/det;
	if(0 <= u && u <= 1){
		*out = GT_Point_add(a_p, GT_Point_scale(a_o, t));
		return 1;
	}
	return 0;
}

size_t GT_Polygon_intersectConvexHalfplane(GT_Point *out, size_t n, const GT_Point points[static n], GT_Halfplane halfplane){
	size_t j = 0, i = 0;
	if(GT_Halfplane_contains(halfplane, points[0])){
		do{
			out[j++] = points[i++];
			if(i == n){
				return n;
			}
		}while(GT_Halfplane_contains(halfplane, points[i]));
		GT_Halfplane_intersectSegment(out + j++, halfplane, points[i - 1], points[i]);
		++i;
		for(; i < n; ++i){
			if(GT_Halfplane_contains(halfplane, points[i])){
				break;
			}
		}
		if(i == n){
			GT_Halfplane_intersectSegment(out + j++, halfplane, points[n - 1], points[0]);
			return j;
		}
		GT_Halfplane_intersectSegment(out + j++, halfplane, points[i - 1], points[i]);
		while(i < n){
			out[j++] = points[i++];
		}
		return j;
	}
	for(i = 1; i < n; ++i){
		if(GT_Halfplane_contains(halfplane, points[i])){
			break;
		}
	}
	if(i == n){
		return 0;
	}
	do{
		out[j++] = points[i++];
	}while(i < n && GT_Halfplane_contains(halfplane, points[i]));
	return j;
}

