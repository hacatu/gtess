#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL2_gfxPrimitives.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <errno.h>
#include <stdio.h>
#include "hash.h"
#include "polygon.h"
#include "polygon_generators.h"
#include "interactive.h"

#define TOOLMENU_LEN 3

int GT_Polygon_append(GT_Polygon *self, GT_Point p){
	if(self->len == self->cap){
		size_t new_cap = self->cap ? self->cap << 1 : 4;
		void *tmp = realloc(self->points, new_cap*sizeof(GT_Point));
		if(!tmp){
			return 0;
		}
		self->points = tmp;
		self->cap = new_cap;
	}
	self->points[self->len++] = p;
	return 1;
}

typedef struct{
	GT_Polygon p;
	size_t bi;
	double A;
} GT_AugPolygon;

struct{
	GT_AugPolygon *ps;
	size_t len, cap;
} polygons;//TODO: should I use a slab allocated linked list?
GT_Interactive ctx;
GT_Hitrect toolMenu_hitrects[TOOLMENU_LEN] = {};

void rotatePointList(size_t n, GT_Point points[static n], off_t k){
	k = (-k)%n;
	if(k < 0){
		k += n;
	}
	size_t orbits = GT_gcd(n, k);
	for(size_t j = 0, i; j < orbits; ++j){
		GT_Point t = points[j];
		for(i = j; i < n - orbits; i += orbits){
			points[i] = points[i + orbits];
		}
		points[i] = t;
	}
}

int canonicalizeConvex(GT_AugPolygon *p){
	printf("Canonicalizing %zu-gon...\n", p->p.len);
	p->A = GT_Polygon_area(p->p.len, p->p.points);
	printf("Signed area is %f\n", p->A);
	if(!GT_Polygon_checkConvex(p->p.len, p->p.points, p->A < 0 ? 1 : 0)){
		printf("The polygon is not convex!\n");
		return 0;
	}
	size_t ai, bi;
	GT_Polygon_diameter(p->p.len, p->p.points, &ai, &bi, p->A < 0 ? 1 : 0);
	printf("The diameter is between points %zu and %zu.\n", ai, bi);
	rotatePointList(p->p.len, p->p.points, -ai);
	if(bi < ai){
		bi = bi + p->p.len - ai;
	}else{
		bi -= ai;
	}
	p->bi = bi;
	printf("Rotated point list so the diameter is between points 0 and %zu.\n", bi);
	if(p->A < 0){
		printf("Flipping orientation...\n");//we could do this in one step but no thanks
		for(size_t i = 1, j = p->p.len - 1; i < j; ++i, --j){
			GT_Point t = p->p.points[i];
			p->p.points[i] = p->p.points[j];
			p->p.points[j] = t;
		}
		//points[i = n - 1 - k] = points[bi = j = 1 + k], k = bi - 1, i = n - bi
		p->bi = p->p.len - p->bi;
		p->A = -p->A;
		printf("Flipped bi is %zu, diagonal between (%f, %f) and (%f, %f)\n", p->bi, p->p.points[0].x, p->p.points[0].y, p->p.points[p->bi].x, p->p.points[p->bi].y);
	}
	return 1;
}

static void init(){
	GT_Interactive_init(&ctx, 1280, 960);
	GT_Interactive_loadTexture(&ctx, "yeet", "/usr/share/icons/HighContrast/scalable/actions/insert-object.svg");
	printf("Textures loaded: {");
	for(GT_Texture *it = hash_next(&ctx.textures, &texture_hashtbl_ft, NULL); it; it = hash_next(&ctx.textures, &texture_hashtbl_ft, it)){
		printf("\"%s\", ", it->name);
	}
	printf("}\n");
}

static void drawToolMenu(int x, int y){
	x = GT_clamp(150, x, ctx.vp.screen_width - 150 - 1);
	y = GT_clamp(150, y, ctx.vp.screen_height - 150 - 1);
	toolMenu_hitrects[0] = (GT_Hitrect){
		.min_x= x - 50, .max_x= x + 50,
		.min_y= y - 150, .max_y= y - 50
	};
	toolMenu_hitrects[1] = (GT_Hitrect){
		.min_x= x - 50, .max_x= x + 50,
		.min_y= y + 50, .max_y= y + 150
	};
	toolMenu_hitrects[2] = (GT_Hitrect){
		.min_x= x - 150, .max_x= x - 50,
		.min_y= y - 100, .max_y= y
	};
	SDL_SetRenderDrawColor(ctx.vp.renderer, 255, 0, 0, 255);
	SDL_Rect rect = GT_Hitrect_toSDL(toolMenu_hitrects[0]);
	SDL_RenderFillRect(ctx.vp.renderer, &rect);
	SDL_SetRenderDrawColor(ctx.vp.renderer, 0, 255, 0, 255);
	rect = GT_Hitrect_toSDL(toolMenu_hitrects[1]);
	SDL_RenderFillRect(ctx.vp.renderer, &rect);
	SDL_SetRenderDrawColor(ctx.vp.renderer, 255, 0, 255, 255);
	rect = GT_Hitrect_toSDL(toolMenu_hitrects[2]);
	SDL_RenderFillRect(ctx.vp.renderer, &rect);
}

static GT_UIState state_drawing(SDL_Event);
static GT_UIState state_toolMenu(SDL_Event);

GT_UIState state_drawing(SDL_Event e){
	GT_UIState next_state = {state_drawing};
	if(e.type == SDL_QUIT){
		GT_Interactive_destroy(&ctx);
		//TODO: we could cleanup polygons, but it's all reachable memory so it doesn't matter
		exit(0);
	}else if(e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT){//TODO: remove extra polygons
		int x, y;
		SDL_GetMouseState(&x, &y);
		GT_Point vptarget = GT_Interactive_fromScreen(&ctx.vp, x, y);
		printf("Adding point at paper coordinates (%f, %f)\n", vptarget.x, vptarget.y);
		GT_Polygon_append(&polygons.ps[0].p, vptarget);//TODO: check for failed alloc
		SDL_SetRenderDrawColor(ctx.vp.renderer, 0x00, 0x00, 0x00, 0xFF);
		SDL_RenderClear(ctx.vp.renderer);
		GT_Interactive_drawPolygonRGBA(&ctx.vp, polygons.ps[0].p, 0xFF, 0x00, 0xFF, 0xFF);
		next_state.dirty = 1;
	}else if(e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_RIGHT){
		int x, y;
		SDL_GetMouseState(&x, &y);
		drawToolMenu(x, y);
		next_state.state_fn = state_toolMenu;
		next_state.dirty = 1;
	}
	return next_state;
}

GT_UIState state_toolMenu(SDL_Event e){
	GT_UIState next_state = {state_toolMenu};
	if(e.type == SDL_QUIT){
		GT_Interactive_destroy(&ctx);
		exit(0);
	}else if(e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT){
		int x, y;
		SDL_GetMouseState(&x, &y);
		GT_Point p = {x, y};
		switch(GT_Hitrect_findTargetIndex(TOOLMENU_LEN, toolMenu_hitrects, p)){
		case 0: {
			polygons.ps[0].A =
				randomConvexNGonWithDiameter(&polygons.ps[0].bi,
					polygons.ps[0].p.len, polygons.ps[0].p.points,
					GT_Point_e1, (GT_Point){-1., 0.}, ctx.rng);
			GT_Interactive_updateContentBounds(&ctx.vp, polygons.ps[0].p.len, polygons.ps[0].p.points);
			SDL_SetRenderDrawColor(ctx.vp.renderer, 0x00, 0x00, 0x00, 0xFF);
			SDL_RenderClear(ctx.vp.renderer);
			GT_Interactive_drawPolygonRGBA(&ctx.vp, polygons.ps[0].p, 0xFF, 0x00, 0xFF, 0xFF);
			next_state.state_fn = state_drawing;
			next_state.dirty = 1;
			break;
		}
		case 1: {
			polygons.ps[0].p.len = 0;
			SDL_SetRenderDrawColor(ctx.vp.renderer, 0x00, 0x00, 0x00, 0xFF);
			SDL_RenderClear(ctx.vp.renderer);
			next_state.state_fn = state_drawing;
			next_state.dirty = 1;
			break;
		}
		case 2: {
			//TODO: ensure the first polygon is convex and ccw.  If it is convex and cw it should just be flipped.
			if(!canonicalizeConvex(&polygons.ps[0])){
				SDL_SetRenderDrawColor(ctx.vp.renderer, 0x00, 0x00, 0x00, 0xFF);
				SDL_RenderClear(ctx.vp.renderer);
				GT_Interactive_drawPolygonRGBA(&ctx.vp, polygons.ps[0].p, 0xFF, 0x00, 0xFF, 0xFF);
				next_state.state_fn = state_drawing;
				next_state.dirty = 1;
			}
			//TODO: how can I pick a good number of points for the intersecting polygon?  3-2n seems reasonable
			size_t n = gsl_rng_uniform_int(ctx.rng, 2*polygons.ps[0].p.len - 3) + 3;
			printf("Generating intersecting %zu-gon...\n", n);
			GT_Point p = randomPointIn(polygons.ps[0].p.len, polygons.ps[0].bi, polygons.ps[0].p.points, polygons.ps[0].A, ctx.rng);
			GT_Point r = randomPointGaussian(p, 2*sqrt(polygons.ps[0].A), ctx.rng);
			printf("First point at (%f, %f)\nSecond point at (%f, %f)\n", p.x, p.y, r.x, r.y);
			if(polygons.ps[1].p.cap < n){
				void *t = realloc(polygons.ps[1].p.points, n*sizeof(GT_Point));//TODO: handle failed alloc
				polygons.ps[1].p.points = t;
				polygons.ps[1].p.cap = n;
			}
			polygons.ps[1].A = randomConvexNGonWithDiameter(&polygons.ps[1].bi, n, polygons.ps[1].p.points, p, r, ctx.rng);
			polygons.ps[1].p.len = n;
			GT_Interactive_updateContentBounds(&ctx.vp, polygons.ps[0].p.len, polygons.ps[0].p.points);
			GT_Interactive_expandContentBounds(&ctx.vp, polygons.ps[1].p.len, polygons.ps[1].p.points);
			SDL_SetRenderDrawColor(ctx.vp.renderer, 0x00, 0x00, 0x00, 0xFF);
			SDL_RenderClear(ctx.vp.renderer);
			GT_Interactive_drawPolygonRGBA(&ctx.vp, polygons.ps[0].p, 0xFF, 0x00, 0xFF, 0xFF);
			GT_Interactive_drawPolygonRGBA(&ctx.vp, polygons.ps[1].p, 0xFF, 0x00, 0xFF, 0xFF);
			//TODO: should I remove the second polygon when done since right non no other state requires it?
			/*
			GT_Point point;
			GT_Polygon dummy = {.points= &point, .len=1, .cap=1};
			gsl_rng *rng_revert = gsl_rng_clone(ctx.rng);
			for(size_t i = 0; i < 200; ++i){
				gsl_rng_memcpy(rng_revert, ctx.rng);
				point = randomPointIn(polygons.ps[0].p.len, polygons.ps[0].bi, polygons.ps[0].p.points, polygons.ps[0].A, ctx.rng);
				GT_Interactive_drawPolygonRGBA(&ctx.vp, dummy, 0x00, 0xFF, 0x00, 0xFF);
				if(!GT_Polygon_contains(polygons.ps[0].p.len, polygons.ps[0].p.points, point, 0)){
					randomPointIn(polygons.ps[0].p.len, polygons.ps[0].bi, polygons.ps[0].p.points, polygons.ps[0].A, rng_revert);
				}
			}
			*/
			next_state.state_fn = state_drawing;
			next_state.dirty = 1;
			break;
		}
		default: {
			break;
		}
		}
	}
	return next_state;
}

int main(){
	polygons.ps = calloc(2, sizeof(*polygons.ps));//TODO: check for failed alloc
	polygons.len = 1;
	polygons.cap = 2;
	init();
	GT_UIState uistate = {.state_fn=state_drawing, .dirty=1};
	for(SDL_Event e;;){
		if(uistate.dirty){
			SDL_RenderPresent(ctx.vp.renderer);
		}
		if(!SDL_WaitEvent(&e)){
			printf("SDL_WaitEvent failed!  SDL Error: %s\n", SDL_GetError());
			exit(EXIT_FAILURE);
		}
		uistate = uistate.state_fn(e);
	}
}

