#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL2_gfxPrimitives.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <errno.h>
#include <stdio.h>
#include "polygon.h"
#include "polygon_generators.h"
#include "interactive.h"

const int SCREEN_WIDTH = 1280;
const int SCREEN_HEIGHT = 960;

SDL_Window *window;
SDL_Renderer *renderer;

void quit(void){
	if(renderer){
		SDL_DestroyRenderer(renderer);
	}
	if(window){
		SDL_DestroyWindow(window);
	}
	IMG_Quit();
	SDL_Quit();
}

void init(void){
	if(SDL_Init(SDL_INIT_VIDEO)){
		printf("SDL could not initialize! SDL Error: %s\n", SDL_GetError());
		exit(EXIT_FAILURE);
	}
	if(!IMG_Init(IMG_INIT_PNG)){
		printf("SDL_image could not initialize! SDL_image ErrorL %s\n", IMG_GetError());
		SDL_Quit();
		exit(EXIT_FAILURE);
	}
	atexit(quit);
	if(!SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "1")){
		printf("Warning: Linear texture filtering not enabled!");
	}
	window = SDL_CreateWindow("GTess Polygon Sandbox", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
	if(!window){
		printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
		exit(EXIT_FAILURE);
	}
	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
	if(renderer == NULL){
		printf("Renderer could not be created! SDL Error: %s\n", SDL_GetError());
		exit(EXIT_FAILURE);
	}
}

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

int main(){
	init();
    gsl_rng_env_setup();
    const gsl_rng_type *rng_T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(rng_T);
	GT_Viewport vp = {.screen_width=SCREEN_WIDTH, .screen_height=SCREEN_HEIGHT, .renderer=renderer, .scale_factor=1.};
	GT_Polygon poly = {};
	int making_shape = 0;
	for(int dirty = 1; 1; dirty = 0){
		for(SDL_Event e; SDL_PollEvent(&e);){
			if(e.type == SDL_QUIT){
				gsl_rng_free(rng);
				free(poly.points);
				exit(0);
			}else if(e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT){
				int x, y;
				SDL_GetMouseState(&x, &y);
				if(!making_shape){
					making_shape = 1;
					poly.len = 0;
				}
				GT_Polygon_append(&poly, GT_Interactive_fromScreen(&vp, x, y));//TODO: check for failed alloc
				SDL_SetRenderDrawColor(renderer, 0x00, 0x00, 0x00, 0xFF);
				SDL_RenderClear(renderer);
				GT_Interactive_drawPolygonRGBA(&vp, poly, 0xFF, 0x00, 0xFF, 0xFF);
				dirty = 1;
			}else if(e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_RIGHT){
				making_shape = 0;
				size_t bi;
				double A = randomConvexNGonWithDiameter(&bi, poly.len, poly.points, GT_Point_e1, (GT_Point){-1., 0.}, rng);
				printf("area by builtin trapezoid decomp: %f\n", A);
				A = GT_Polygon_area(poly.len, poly.points);
				printf("area by shoestring method: %f\n", A);
				GT_Interactive_updateContentBounds(&vp, poly.len, poly.points);
				SDL_SetRenderDrawColor(renderer, 0x00, 0x00, 0x00, 0xFF);
				SDL_RenderClear(renderer);
				GT_Interactive_drawPolygonRGBA(&vp, poly, 0xFF, 0x00, 0xFF, 0xFF);
				dirty = 1;
			}
		}
		if(dirty){
			SDL_RenderPresent(renderer);
		}
	}
}

