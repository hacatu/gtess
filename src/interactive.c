#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL2_gfxPrimitives.h>
#include <errno.h>
#include <stdio.h>
#include "polygon.h"
#include "interactive.h"

void GT_Interactive_updateContentBounds(GT_Viewport *self, size_t len, const GT_Point points[static len]){
	if(!len){
		return;
	}
	self->min_x = self->max_x = points[0].x;
	self->min_y = self->max_y = points[0].y;
	for(size_t i = 1; i < len; ++i){
		if(points[i].x < self->min_x){
			self->min_x = points[i].x;
		}else if(points[i].x > self->max_x){
			self->max_x = points[i].x;
		}
		if(points[i].y < self->min_y){
			self->min_y = points[i].y;
		}else if(points[i].y > self->max_y){
			self->max_y = points[i].y;
		}
	}
	GT_Interactive_updateScaleFactor(self);
}

void GT_Interactive_updateScaleFactor(GT_Viewport *self){
	const double content_width = self->max_x - self->min_x;
	const double content_height = self->max_y - self->min_y;
	if(content_width < GT_EPSILON){
		if(content_height < GT_EPSILON){
			self->scale_factor = 1.;
		}else{
			self->scale_factor = self->screen_height/content_height;
		}
	}else if(content_height < GT_EPSILON){
		self->scale_factor = self->screen_width/content_width;
	}else{
		self->scale_factor = fmin(self->screen_width/content_width, self->screen_height/content_height);
	}
}

GT_Point GT_Interactive_toScreen(const GT_Viewport *self, GT_Point p){
	return (GT_Point){
		self->scale_factor*(p.x - self->min_x),
		self->scale_factor*(self->max_y - p.y)
	};
}

GT_Point GT_Interactive_fromScreen(const GT_Viewport *self, int x, int y){
	return (GT_Point){
		x/self->scale_factor + self->min_x,
		self->max_y - y/self->scale_factor
	};
}

void GT_Interactive_drawPolygonRGBA(GT_Viewport *vp, GT_Polygon poly, Uint8 r, Uint8 g, Uint8 b, Uint8 a){
		for(size_t i = 0; i < poly.len; ++i){
			GT_Point p = GT_Interactive_toScreen(vp, poly.points[i]);
			aacircleRGBA(vp->renderer, p.x, p.y, 5, r, g, b, a);
		}
		for(size_t i = 0; i + 1 < poly.len; ++i){
			GT_Point p = GT_Interactive_toScreen(vp, poly.points[i]);
			GT_Point q = GT_Interactive_toScreen(vp, poly.points[i + 1]);
			aalineRGBA(vp->renderer, p.x, p.y, q.x, q.y, r, g, b, a);
		}
		GT_Point p = GT_Interactive_toScreen(vp, poly.points[poly.len - 1]);
		GT_Point q = GT_Interactive_toScreen(vp, poly.points[0]);
		aalineRGBA(vp->renderer, p.x, p.y, q.x, q.y, r, g, b, a);
}

