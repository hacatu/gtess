#ifndef GTESS_INTERACTIVE_H
#define GTESS_INTERACTIVE_H

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL2_gfxPrimitives.h>
#include <errno.h>
#include <stdio.h>
#include "hash.h"
#include "polygon.h"

typedef struct{
	double min_x, max_x, min_y, max_y;
	void *ent;
} GT_Hitrect;

typedef struct{
	int screen_width, screen_height;
	double min_x, max_x, min_y, max_y;
	SDL_Renderer *renderer;
	double scale_factor;
} GT_Viewport;

typedef struct{
	const char *name;
	SDL_Texture *tex;
} GT_Texture;

typedef struct GT_Interactive GT_Interactive;
struct GT_Interactive{
	SDL_Window *window;
	GT_Viewport vp;
	gsl_rng *rng;
	hashtbl_t textures;
};

typedef struct GT_UIState GT_UIState;
struct GT_UIState{
	GT_UIState (*state_fn)();
	int dirty;
};

void GT_Interactive_init(GT_Interactive *self, int screen_width, int screen_height);
void GT_Interactive_destroy(GT_Interactive *self);
int GT_Interactive_loadTexture(GT_Interactive *self, const char *name, const char *path);
SDL_Texture *GT_Interactive_getTexture(GT_Interactive *self, const char *name);

int GT_Hitrect_contains(GT_Hitrect hitrect, GT_Point p);
//returns len on not found
size_t GT_Hitrect_findTargetIndex(size_t len, const GT_Hitrect hitrects[static len], GT_Point p);
void *GT_Hitrect_findTargetEnt(size_t len, const GT_Hitrect hitrects[static len], GT_Point p);

void GT_Interactive_updateContentBounds(GT_Viewport *self, size_t len, const GT_Point points[static len]);
void GT_Interactive_updateScaleFactor(GT_Viewport *self);

GT_Point GT_Interactive_toScreen(const GT_Viewport *self, GT_Point p);
GT_Point GT_Interactive_fromScreen(const GT_Viewport *self, int x, int y);

void GT_Interactive_drawPolygonRGBA(GT_Viewport *vp, GT_Polygon poly, Uint8 r, Uint8 g, Uint8 b, Uint8 a);

extern hashtbl_ft texture_hashtbl_ft;

#endif

