#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "point.h"
#include "region.h"
#include "paper_initializers.h"

GT_TwistNode *GT_makeExterior_square(double side){
    GT_TwistNode *exterior_node = malloc(offsetof(GT_TwistNode, data_buf) + 4*(2*sizeof(GT_Point) + sizeof(GT_PleatProfile) + sizeof(GT_TaggedNodePtr)));
    if(!exterior_node){
        return NULL;
    }
    GT_FlatNode *interior_node = malloc(offsetof(GT_FlatNode, data_buf) + 4*(sizeof(GT_Point) + sizeof(GT_TaggedNodePtr)));
    exterior_node->n_sides = 4;
    exterior_node->scratch_data = 0;
    exterior_node->is_external = 1;
    GT_Point *exterior_vertices = GT_TwistRegion_getVertices(exterior_node);
    exterior_vertices[0] = (GT_Point){0, 0};
    exterior_vertices[1] = (GT_Point){0, side};
    exterior_vertices[2] = (GT_Point){side, side};
    exterior_vertices[3] = (GT_Point){side, 0};
    GT_Point *exterior_directions = GT_TwistRegion_getDirections(exterior_node);
    exterior_directions[0] = (GT_Point){1, 0};
    exterior_directions[1] = (GT_Point){0, -1};
    exterior_directions[2] = (GT_Point){-1, 0};
    exterior_directions[3] = (GT_Point){0, 1};
    memset(GT_TwistRegion_getPleatProfiles(exterior_node), 0, 4*sizeof(GT_PleatProfile));
    GT_TaggedNodePtr *exterior_adjs = GT_TwistNode_getAdjPleats(exterior_node);
    for(size_t i = 0; i < 4; ++i){
        exterior_adjs[i] = (GT_TaggedNodePtr){.tag=GT_FlatTag, .flat=interior_node, .back_i=3 - i};
    }
    interior_node->n_sides = 4;
    interior_node->scratch_data = 0;
    GT_Point *interior_vertices = GT_FlatRegion_getVertices(interior_node);
    for(size_t i = 0; i < 4; ++i){
        interior_vertices[i] = exterior_vertices[(4 - i)%4];
    }
    GT_TaggedNodePtr *interior_adjs = GT_FlatNode_getAdjPleats(interior_node);
    for(size_t i = 0; i < 4; ++i){
        interior_adjs[i] = (GT_TaggedNodePtr){.tag=GT_TwistTag, .twist=exterior_node, .back_i=3 - i};
    }
    return exterior_node;
}

