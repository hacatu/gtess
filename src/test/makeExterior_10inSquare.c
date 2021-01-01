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

static int checkTwistRegion(const GT_TwistRegion *self, int print_err){
    if(!GT_Polygon_checkConvex(self->n_sides, self->vertices, self->is_external)){
        DEBUG_PRINTF(print_err, "Twist node %p is't convex\n", self);
        return 0;
    }else if(!GT_Polygon_checkOrigin_strided(self->n_sides, sizeof(GT_TwistNeighbor), &self->neighbors[0].direction, self->is_external)){
        DEBUG_PRINTF(print_err, "Twist node %p doesn't have counterclockwise normals\n", self);
        return 0;
    }
    for(size_t i = 0; i < self->n_sides; ++i){
        GT_TaggedRegionPtr region = self->neighbors[i].region;
        if(self->is_external && region.tag == GT_FlatTag){
            if(!GT_checkTwistFlatAdj(self, i, region.flat, region.back_i)){
                DEBUG_PRINTF(print_err, "Twist node %p isn't adjacent to flat node %p\n", self, region.flat);
                return 0;
            }
        }else if(!GT_checkTwistPleatAdj(self, i, region.pleat, region.back_i)){
            DEBUG_PRINTF(print_err, "Twist node %p isn't adjacent to pleat node %p\n", self, region.pleat);
            return 0;
        }
    }
    return 1;
}

static int checkPleatRegion(const GT_PleatRegion *self, int print_err){
    if(!GT_Polygon_checkConvex(self->n_sides, self->vertices, 0)){
        DEBUG_PRINTF(print_err, "Pleat node %p isn't convex\n", self);
        return 0;
    }else if(!GT_PleatRegion_checkProfiles(self)){
        DEBUG_PRINTF(print_err, "Pleat node %p doesn't have consistent pleat profiles\n", self);
        return 0;
    }
    for(size_t i = 0; i < self->n_sides; ++i){
        GT_TaggedRegionPtr region = self->neighbors[i].region;
        if(region.tag == GT_TwistTag){
            if(!GT_checkTwistPleatAdj(region.twist, region.back_i, self, i)){
                DEBUG_PRINTF(print_err, "Pleat node %p isn't adjacent to twist node %p\n", self, region.twist);
                return 0;
            }
        }else if(!GT_checkFlatPleatAdj(region.flat, region.back_i, self, i)){
            DEBUG_PRINTF(print_err, "Pleat node %p isn't adjacent to flat node %p\n", self, region.twist);
            return 0;
        }
    }
    return 1;
}

static int checkFlatRegion(const GT_FlatRegion *self, int print_err){
    if(!GT_Polygon_checkConvex(self->n_sides, self->vertices, 0)){
        DEBUG_PRINTF(print_err, "Flat node %p isn't convex\n", self);
        return 0;
    }
    for(size_t i = 0; i < self->n_sides; ++i){
        GT_TaggedRegionPtr region = self->neighbors[i].region;
        if(region.tag == GT_PleatTag){
            if(!GT_checkFlatPleatAdj(self, i, region.pleat, region.back_i)){
                DEBUG_PRINTF(print_err, "Flat node %p isn't adjacent to pleat node %p\n", self, region.pleat);
                return 0;
            }
        }else if(!GT_checkTwistFlatAdj(region.twist, region.back_i, self, i)){
            DEBUG_PRINTF(print_err, "Flat node %p isn't adjacent to twist node %p\n", self, region.twist);
            return 0;
        }
    }
    return 1;
}

static void checkRegion_cb(GT_TaggedRegionPtr region, void *data){
    switch(region.tag){
        case GT_TwistTag:
            *(int*)data &= checkTwistRegion(region.twist, 1);
            break;
        case GT_PleatTag:
            *(int*)data &= checkPleatRegion(region.pleat, 1);
            break;
        case GT_FlatTag:
            *(int*)data &= checkFlatRegion(region.flat, 1);
            break;
    }

}

static int checkInternals(const GT_TwistRegion *exterior){
    const static GT_Point expected_exterior_vertices[] = {{0, 0}, {0, 10}, {10, 10}, {10, 0}};
    const static GT_Point expected_exterior_directions[] = {{1, 0}, {0, -1}, {-1, 0}, {0, 1}};
    const static GT_Point expected_interior_vertices[] = {{0, 0}, {10, 0}, {10, 10}, {0, 10}};

    if(exterior->n_sides != 4 || !exterior->is_external){
        return 0;
    }

    for(size_t i = 0; i < exterior->n_sides; ++i){
        if(exterior->vertices[i].x != expected_exterior_vertices[i].x ||
                exterior->vertices[i].y != expected_exterior_vertices[i].y ||
                exterior->neighbors[i].direction.x != expected_exterior_directions[i].x ||
                exterior->neighbors[i].direction.y != expected_exterior_directions[i].y ||
                exterior->neighbors[i].profile.n_pleats || exterior->neighbors[i].region.tag != GT_FlatTag){
            return 0;
        }
    }
    const GT_FlatRegion *interior = exterior->neighbors[0].region.flat;
    if(interior->n_sides != exterior->n_sides || interior->scratch_data != exterior->scratch_data){
        return 0;
    }
    for(size_t i = 0; i < exterior->n_sides; ++i){
        if(interior->vertices[i].x != expected_interior_vertices[i].x ||
                interior->vertices[i].y != expected_interior_vertices[i].y ||
                interior->neighbors[i].region.tag != GT_TwistTag){
            return 0;
        }
    }
    for(size_t i = 0; i < exterior->n_sides; ++i){
        if(exterior->neighbors[i].region.flat != interior || exterior->neighbors[i].region.back_i != 3 - i){
            return 0;
        }else if(interior->neighbors[i].region.twist != exterior || interior->neighbors[i].region.back_i != 3 - i){
            return 0;
        }
    }
    return 1;
}

int main(){
    GT_TwistRegion *exterior = GT_makeExterior_square(10);
    if(!exterior){
        printf("\e[1;31mERROR: Could not create external node.\e[0m\n");
        exit(EXIT_FAILURE);
    }
    if(!checkInternals(exterior)){
        printf("\e[1;31mERROR: Initial configuration does not have the expected internal layout.\e[0m\n");
        exit(EXIT_FAILURE);
    }
    printf("exterior: %p\n", exterior);
    printf("interior: %p\n", exterior->neighbors[0].region.flat);
    int passed = 1;
    GT_forEachPaperRegion(exterior, &passed, checkRegion_cb);
    if(!passed){
        printf("\e[1;31mERROR: Initial configuration failed consistency checks.\e[0m\n");
        exit(EXIT_FAILURE);
    }
    GT_deletePaper(exterior);
}

