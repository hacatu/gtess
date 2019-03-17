#ifndef GTESS_INTERACTIVE_H
#define GTESS_INTERACTIVE_H

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL2_gfxPrimitives.h>
#include <errno.h>
#include <stdio.h>
#include "polygon.h"

typedef struct{
	int screen_width, screen_height;
	double min_x, max_x, min_y, max_y;
	SDL_Renderer *renderer;
	double scale_factor;
} GT_Viewport;

void GT_Interactive_updateContentBounds(GT_Viewport *self, size_t len, const GT_Point points[static len]);
void GT_Interactive_updateScaleFactor(GT_Viewport *self);

GT_Point GT_Interactive_toScreen(const GT_Viewport *self, GT_Point p);
GT_Point GT_Interactive_fromScreen(const GT_Viewport *self, int x, int y);

void GT_Interactive_drawPolygonRGBA(GT_Viewport *vp, GT_Polygon poly, Uint8 r, Uint8 g, Uint8 b, Uint8 a);

#endif

