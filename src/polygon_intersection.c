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

static inline size_t gcd(size_t a, size_t b){
	if(a < b){
		b %= a;
	}
	while(1){
		if(!b){
			return a;
		}
		a %= b;
		if(!a){
			return b;
		}
		b %= a;
	}
}

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
		size_t orbits = gcd(self->cap, self->cap - self->a);
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
	size_t step = is_top ? 1 : polygon->n - 1;
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
		}//we don't need to check for overlap anymore because overlap of an intersection of open segments implies the end of one of the polygons
	}
	return 1;
}

static inline int GT_SLE_scan_stateless(GT_SLE_State *state){
	GT_SLE_Event e;
	GT_EventHeap_pop(&state->events, &e);
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
	if(min_overlap <= max_overlap){
		if(fabs(e.p.y - min_overlap) < GT_EPSILON){
			GT_Point_Queue_pusha(&state->points, e.p);
		}else{
			GT_Point_Queue_pushb(&state->points, e.p);
		}
	}else if(state->points.len){
		return 0;
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
	return GT_Polygon_area(polygon->n, polygon->points) > GT_EPSILON;//if it is smaller the polygon is either a line or listed clockwise
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

size_t GT_Polygon_intersectConvex(GT_Point *points_c, size_t n_a, const GT_Point points_a[static n_a], size_t n_b, const GT_Point points_b[static n_b]){
	GT_SLE_Event _event_buf[8];
	GT_SLE_State state = {
		.polygons[0]={.points=points_a, .n=n_a},
		.polygons[1]={.points=points_b, .n=n_b},
		.points={.buf=points_c, .cap= n_a + n_b},
		.events={.buf=_event_buf, .cap=8},
		.state_fn=GT_SLE_scan_stateless,
	};
	int a_line = !GT_SLE_prepare_polygon(state.polygons);
	int b_line = !GT_SLE_prepare_polygon(state.polygons + 1);
	if(a_line || b_line){
		//handle simple case
	}
	if(points_b[state.polygons[1].seg_ends[1][0]].x < points_a[state.polygons[0].seg_ends[1][0]].x){
		GT_SLE_Polygon t = state.polygons[0];
		state.polygons[0] = state.polygons[1];
		state.polygons[1] = t;
	}
	if(!GT_SLE_scan_to_start(&state, 0)){
		return 0;
	}
	GT_SLE_scan_to_start(&state, 1);
	while(state.state_fn(&state));
	GT_Point_Queue_canonicalize(&state.points);
	return state.points.len;
}

