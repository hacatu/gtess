#define _GNU_SOURCE
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "point.h"
#include "polygon.h"
#include "region.h"
#include "paper_initializers.h"

#define DEBUG_PRINTF(do_print, printf_args...) (do_print ? printf(printf_args) : 0)

static int checkTwistNode(const GT_TwistNode *self, int print_err){
    const GT_Point *vertices = GT_TwistRegion_getVertices((GT_TwistRegion*)self);
    const GT_Point *directions = GT_TwistRegion_getDirections((GT_TwistRegion*)self);
    if(!GT_Polygon_checkConvex(self->n_sides, vertices, self->is_external)){
        DEBUG_PRINTF(print_err, "Twist node %p is't convex\n", self);
        return 0;
    }else if(!GT_Polygon_checkOrigin(self->n_sides, directions, self->is_external)){
        DEBUG_PRINTF(print_err, "Twist node %p doesn't have counterclockwise normals\n", self);
        return 0;
    }
    const GT_TaggedNodePtr *adjs = GT_TwistNode_getAdjPleats((GT_TwistRegion*)self);
    for(size_t i = 0; i < self->n_sides; ++i){
        GT_TaggedNodePtr tnp = adjs[i];
        if(self->is_external && tnp.tag == GT_FlatTag){
            if(!GT_checkTwistFlatAdj(self, i, tnp.flat, tnp.back_i)){
                DEBUG_PRINTF(print_err, "Twist node %p isn't adjacent to flat node %p\n", self, tnp.flat);
                return 0;
            }
        }else if(!GT_checkTwistPleatAdj(self, i, tnp.pleat, tnp.back_i)){
            DEBUG_PRINTF(print_err, "Twist node %p isn't adjacent to pleat node %p\n", self, tnp.pleat);
            return 0;
        }
    }
    return 1;
}

static int checkPleatNode(const GT_PleatNode *self, int print_err){
    const GT_Point *vertices = GT_PleatRegion_getVertices((GT_PleatRegion*)self);
    if(!GT_Polygon_checkConvex(self->n_sides, vertices, 0)){
        DEBUG_PRINTF(print_err, "Pleat node %p isn't convex\n", self);
        return 0;
    }else if(!GT_PleatRegion_checkProfiles(self)){
        DEBUG_PRINTF(print_err, "Pleat node %p doesn't have consistent pleat profiles\n", self);
        return 0;
    }
    const GT_TaggedNodePtr *adjs = GT_PleatNode_getAdjNodes((GT_PleatRegion*)self);
    for(size_t i = 0; i < self->n_sides; ++i){
        GT_TaggedNodePtr tnp = adjs[i];
        if(tnp.tag == GT_TwistTag){
            if(!GT_checkTwistPleatAdj(tnp.twist, tnp.back_i, self, i)){
                DEBUG_PRINTF(print_err, "Pleat node %p isn't adjacent to twist node %p\n", self, tnp.twist);
                return 0;
            }
        }else if(!GT_checkFlatPleatAdj(tnp.flat, tnp.back_i, self, i)){
            DEBUG_PRINTF(print_err, "Pleat node %p isn't adjacent to flat node %p\n", self, tnp.twist);
            return 0;
        }
    }
    return 1;
}

static int checkFlatNode(const GT_FlatNode *self, int print_err){
    const GT_Point *vertices = GT_FlatRegion_getVertices((GT_FlatRegion*)self);
    if(!GT_Polygon_checkConvex(self->n_sides, vertices, 0)){
        DEBUG_PRINTF(print_err, "Flat node %p isn't convex\n", self);
        return 0;
    }
    const GT_TaggedNodePtr *adjs = GT_FlatNode_getAdjPleats((GT_FlatNode*)self);
    for(size_t i = 0; i < self->n_sides; ++i){
        GT_TaggedNodePtr tnp = adjs[i];
        if(tnp.tag == GT_PleatTag){
            if(!GT_checkFlatPleatAdj(self, i, tnp.pleat, tnp.back_i)){
                DEBUG_PRINTF(print_err, "Flat node %p isn't adjacent to pleat node %p\n", self, tnp.pleat);
                return 0;
            }
        }else if(!GT_checkTwistFlatAdj(tnp.twist, tnp.back_i, self, i)){
            DEBUG_PRINTF(print_err, "Flat node %p isn't adjacent to twist node %p\n", self, tnp.twist);
            return 0;
        }
    }
    return 1;
}

static void checkNode_cb(GT_TaggedNodePtr tnp, void *data){
    switch(tnp.tag){
        case GT_TwistTag:
            *(int*)data &= checkTwistNode(tnp.twist, 1);
            break;
        case GT_PleatTag:
            *(int*)data &= checkPleatNode(tnp.pleat, 1);
            break;
        case GT_FlatTag:
            *(int*)data &= checkFlatNode(tnp.flat, 1);
            break;
    }

}

static int checkInternals(const GT_TwistNode *exterior_node){
    const static GT_Point expected_exterior_vertices[] = {{0, 0}, {0, 10}, {10, 10}, {10, 0}};
    const static GT_Point expected_exterior_directions[] = {{1, 0}, {0, -1}, {-1, 0}, {0, 1}};
    const static GT_Point expected_interior_vertices[] = {{0, 0}, {10, 0}, {10, 10}, {0, 10}};

    if(exterior_node->n_sides != 4 || !exterior_node->is_external){
        return 0;
    }

    const GT_Point *exterior_vertices = GT_TwistRegion_getVertices((GT_TwistRegion*)exterior_node);
    const GT_Point *exterior_directions = GT_TwistRegion_getDirections((GT_TwistRegion*)exterior_node);
    const GT_PleatProfile *exterior_pleats = GT_TwistRegion_getPleatProfiles((GT_TwistRegion*)exterior_node);
    const GT_TaggedNodePtr *exterior_adjs = GT_TwistNode_getAdjPleats((GT_TwistNode*)exterior_node);

    for(size_t i = 0; i < exterior_node->n_sides; ++i){
        if(exterior_vertices[i].x != expected_exterior_vertices[i].x ||
                exterior_vertices[i].y != expected_exterior_vertices[i].y ||
                exterior_directions[i].x != expected_exterior_directions[i].x ||
                exterior_directions[i].y != expected_exterior_directions[i].y ||
                exterior_pleats[i].n_pleats || exterior_adjs[i].tag != GT_FlatTag){
            return 0;
        }
    }
    const GT_FlatNode *interior_node = exterior_adjs[0].flat;
    if(interior_node->n_sides != exterior_node->n_sides || interior_node->scratch_data != exterior_node->scratch_data){
        return 0;
    }
    const GT_Point *interior_vertices = GT_FlatRegion_getVertices((GT_FlatRegion*)interior_node);
    const GT_TaggedNodePtr *interior_adjs = GT_FlatNode_getAdjPleats((GT_FlatNode*)interior_node);
    for(size_t i = 0; i < exterior_node->n_sides; ++i){
        if(interior_vertices[i].x != expected_interior_vertices[i].x ||
                interior_vertices[i].y != expected_interior_vertices[i].y ||
                interior_adjs[i].tag != GT_TwistTag){
            return 0;
        }
    }
    for(size_t i = 0; i < exterior_node->n_sides; ++i){
        if(exterior_adjs[i].flat != interior_node || exterior_adjs[i].back_i != 3 - i){
            return 0;
        }else if(interior_adjs[i].twist != exterior_node || interior_adjs[i].back_i != 3 - i){
            return 0;
        }
    }
    return 1;
}

int main(){
    GT_TwistNode *exterior_node = GT_makeExterior_square(10);
    if(!exterior_node){
        printf("\e[1;31mERROR: Could not create external node.\e[0m\n");
        exit(EXIT_FAILURE);
    }
    if(!checkInternals(exterior_node)){
        printf("\e[1;31mERROR: Initial configuration does not have the expected internal layout.\e[0m\n");
        exit(EXIT_FAILURE);
    }
    printf("exterior_node: %p\n", exterior_node);
    printf("interior_node: %p\n", GT_TwistNode_getAdjPleats(exterior_node)->flat);
    int passed = 1;
    GT_forEachPaperNode(exterior_node, &passed, checkNode_cb);
    if(!passed){
        printf("\e[1;31mERROR: Initial configuration failed consistency checks.\e[0m\n");
        exit(EXIT_FAILURE);
    }
    GT_deletePaper(exterior_node);
}

