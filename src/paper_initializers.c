#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "paper_initializers.h"

GT_TwistRegion *GT_makeExterior_square(double side){
    GT_TwistRegion *exterior = malloc(1*sizeof(GT_TwistRegion));
    GT_FlatRegion *interior = malloc(1*sizeof(GT_FlatRegion));
    if(!exterior || !interior){
		free(interior);
		free(exterior);
        return NULL;
    }
    exterior->vertices = malloc(4*sizeof(GT_Point));
    interior->vertices = malloc(4*sizeof(GT_Point));
    exterior->neighbors = malloc(4*sizeof(GT_TwistNeighbor));
    interior->neighbors = malloc(4*sizeof(GT_FlatNeighbor));
    if(!exterior->vertices || !interior->vertices || !exterior->neighbors || !interior->neighbors){
		free(interior->neighbors);
		free(exterior->neighbors);
		free(interior->vertices);
		free(exterior->vertices);
		free(interior);
		free(exterior);
		return NULL;
	}
    exterior->n_sides = 4;
    exterior->scratch_data = 0;
    exterior->is_external = 1;
    for(size_t i = 0; i < 4; ++i){
		exterior->vertices[i] = (GT_Point){i&2 ? side : 0, (i + 1)&2 ? side : 0};
		exterior->neighbors[i] = (GT_TwistNeighbor){
			.direction= GT_Point_dir[(4 - i)%4],
			.region= {.tag=GT_FlatTag, .flat=interior, .back_i= 3 - i}
		};
	}
    interior->n_sides = 4;
    interior->scratch_data = 0;
    for(size_t i = 0; i < 4; ++i){
		interior->vertices[i] = exterior->vertices[(4 - i)%4];
		interior->neighbors[i] = (GT_FlatNeighbor){
			.region= {.tag=GT_TwistTag, .twist=exterior, .back_i= 3 - i}
		};
    }
    return exterior;
}

