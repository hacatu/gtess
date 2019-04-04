#define _GNU_SOURCE
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL2_gfxPrimitives.h>
#include <errno.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "hash.h"
#include "polygon.h"
#include "interactive.h"

const int TEXTURE_HASHTBL_RESERVE = 32;

void GT_Interactive_init(GT_Interactive *self, int screen_width, int screen_height){
	self->vp = (GT_Viewport){.screen_width=screen_width, .screen_height=screen_height, .scale_factor=1.};
	if(SDL_Init(SDL_INIT_VIDEO)){
		printf("SDL could not initialize! SDL Error: %s\n", SDL_GetError());
		exit(EXIT_FAILURE);
	}
	if(!IMG_Init(IMG_INIT_PNG)){
		printf("SDL_image could not initialize! SDL_image Error: %s\n", IMG_GetError());
		SDL_Quit();
		exit(EXIT_FAILURE);
	}
	if(!SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "1")){
		printf("Warning: Linear texture filtering not enabled!");
	}
	self->window = SDL_CreateWindow("GTess Polygon Sandbox", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, self->vp.screen_width, self->vp.screen_height, SDL_WINDOW_SHOWN);
	if(!self->window){
		printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
		exit(EXIT_FAILURE);
	}
	self->vp.renderer = SDL_CreateRenderer(self->window, -1, SDL_RENDERER_ACCELERATED);
	if(self->vp.renderer == NULL){
		printf("Renderer could not be created! SDL Error: %s\n", SDL_GetError());
		exit(EXIT_FAILURE);
	}
    gsl_rng_env_setup();
    const gsl_rng_type *rng_T = gsl_rng_default;
    self->rng = gsl_rng_alloc(rng_T);
    if(!hash_init(&self->textures, &texture_hashtbl_ft, TEXTURE_HASHTBL_RESERVE)){
		printf("Texture hash table could not be created!\n");
		exit(EXIT_FAILURE);
	}
}

void GT_Interactive_destroy(GT_Interactive *self){
	SDL_DestroyRenderer(self->vp.renderer);
	SDL_DestroyWindow(self->window);
	IMG_Quit();
	SDL_Quit();
	gsl_rng_free(self->rng);
	hash_destroy(&self->textures, &texture_hashtbl_ft);
}

int GT_Interactive_loadTexture(GT_Interactive *self, const char *name, const char *path){
	SDL_Surface *img = IMG_Load(path);
	if(!img){
		printf("Image file \"%s\" could not be loaded.  IMG Error: %s\n", path, IMG_GetError());
		return 1;
	}
	SDL_Texture *tex = SDL_CreateTextureFromSurface(self->vp.renderer, img);
	SDL_FreeSurface(img);
	if(!tex){
		printf("Image \"%s\" could not be converted to a texture.  SDL Error: %s\n", name, SDL_GetError());
		return 1;
	}
	int ret;
	hash_insert(&self->textures, &texture_hashtbl_ft, &(GT_Texture){.name=strdup(name), .tex=tex}, &ret);
	return ret;
}

SDL_Texture *GT_Interactive_getTexture(GT_Interactive *self, const char *name){
	GT_Texture *slot = hash_get(&self->textures, &texture_hashtbl_ft, &(GT_Texture){.name=name});
	if(slot){
		return slot->tex;
	}
	return NULL;
}

int GT_Hitrect_contains(GT_Hitrect hitrect, GT_Point p){
	return (p.x - hitrect.min_x <= GT_EPSILON) &&
		(hitrect.max_x - p.x <= GT_EPSILON) &&
		(p.y - hitrect.min_y <= GT_EPSILON) &&
		(hitrect.max_y - p.y <= GT_EPSILON);
}

size_t GT_Hitrect_findTargetIndex(size_t len, const GT_Hitrect hitrects[static len], GT_Point p){
	for(size_t i = 0; i < len; ++i){
		if(GT_Hitrect_contains(hitrects[i], p)){
			return i;
		}
	}
	return len;
}

void *GT_Hitrect_findTargetEnt(size_t len, const GT_Hitrect hitrects[static len], GT_Point p){
	for(size_t i = 0; i < len; ++i){
		if(GT_Hitrect_contains(hitrects[i], p)){
			return hitrects[i].ent;
		}
	}
	return NULL;
}

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

uint64_t texture_hash_fn(const void *e){
	uint64_t h = 0;
	for(const char *s = ((const GT_Texture*)e)->name; *s; ++s){
		h = (h << 5) - h + *s;
	}
	return h;
}

int texture_cmp_fn(const void *a, const void *b){
	return strcmp(((const GT_Texture*)a)->name, ((const GT_Texture*)b)->name);
}

void texture_del_fn(void *a){
	free((char*)((GT_Texture*)a)->name);//discard const
	SDL_DestroyTexture(((GT_Texture*)a)->tex);
}

hashtbl_ft texture_hashtbl_ft = {
	.size= sizeof(GT_Texture),
	.hash= texture_hash_fn,
	.cmp= texture_cmp_fn,
	.add= 0,
	.del= texture_del_fn,
	.load_factor= 0.77
};

