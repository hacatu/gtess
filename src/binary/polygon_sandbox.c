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

#define TOOLMENU_LEN 2

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

GT_Polygon poly;
int making_shape;
GT_Interactive ctx;
GT_Hitrect toolMenu_hitrects[TOOLMENU_LEN] = {};

static void init(){
	GT_Interactive_init(&ctx, 1280, 960);
	GT_Interactive_loadTexture(&ctx, "yeet", "/usr/share/icons/HighContrast/scalable/actions/insert-object.svg");
	printf("Textures loaded: {");
	for(GT_Texture *it = hash_next(&ctx.textures, &texture_hashtbl_ft, NULL); it; it = hash_next(&ctx.textures, &texture_hashtbl_ft, it)){
		printf("\"%s\", ", it->name);
	}
	printf("}\n");
}

static GT_UIState state_drawing(SDL_Event);
static GT_UIState state_toolMenu(SDL_Event);

GT_UIState state_drawing(SDL_Event e){
	GT_UIState next_state = {state_drawing};
	if(e.type == SDL_QUIT){
		GT_Interactive_destroy(&ctx);
		free(poly.points);
		exit(0);
	}else if(e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT){
		int x, y;
		SDL_GetMouseState(&x, &y);
		if(!making_shape){
			making_shape = 1;
			poly.len = 0;
		}
		GT_Polygon_append(&poly, GT_Interactive_fromScreen(&ctx.vp, x, y));//TODO: check for failed alloc
		SDL_SetRenderDrawColor(ctx.vp.renderer, 0x00, 0x00, 0x00, 0xFF);
		SDL_RenderClear(ctx.vp.renderer);
		GT_Interactive_drawPolygonRGBA(&ctx.vp, poly, 0xFF, 0x00, 0xFF, 0xFF);
		next_state.dirty = 1;
	}else if(e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_RIGHT){//TODO: actully what we should do here is pop open the tool menu
		making_shape = 0;
		next_state.state_fn = state_toolMenu;
		size_t bi;
		double A = randomConvexNGonWithDiameter(&bi, poly.len, poly.points, GT_Point_e1, (GT_Point){-1., 0.}, ctx.rng);
		GT_Interactive_updateContentBounds(&ctx.vp, poly.len, poly.points);
		SDL_SetRenderDrawColor(ctx.vp.renderer, 0x00, 0x00, 0x00, 0xFF);
		SDL_RenderClear(ctx.vp.renderer);
		GT_Interactive_drawPolygonRGBA(&ctx.vp, poly, 0xFF, 0x00, 0xFF, 0xFF);
		next_state.dirty = 1;
	}
	return next_state;
}

GT_UIState state_toolMenu(SDL_Event e){
	GT_UIState next_state = {state_toolMenu};
	//TODO: actually check for lclicks on the tools using the toolMenu_hitrects list that state_drawing should populate on rclick
	return next_state;
}

int main(){
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

