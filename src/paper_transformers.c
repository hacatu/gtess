#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "paper_transformers.h"
#include "polygon.h"

static GT_TaggedRegionPtr findPaperRegion_rec(GT_TwistRegion *currTwist, size_t cwi, size_t ccwi, int scratch, GT_Point a);
static GT_TaggedRegionPtr findPaperRegion_base(GT_PleatRegion *pleat, size_t i, GT_Point a);

GT_TaggedRegionPtr GT_findPaperRegion(GT_TwistRegion *external, GT_Point a){
	GT_TaggedRegionPtr reg;
	int scratch = external->scratch_data;
	external->scratch_data = !scratch;
	size_t i;
	for(i = 0; i < external->n_sides; ++i){
		reg = external->neighbors[i].region;
		if(reg.tag == GT_PleatTag){
			break;
		}
	}
	if(i == external->n_sides){//paper is just one flat region
		return reg;
	}
	GT_PleatRegion *pleat = reg.pleat;
	if(GT_Polygon_contains(pleat->n_sides, pleat->vertices, a, 0)){
		return reg;
	}
	size_t external_i = reg.back_i;
	for(size_t i = 0; i < pleat->n_sides; ++i){//technically we are guaranteed to find an internal twist node adjacent to pleat but check bounds anyway
		if(i == external_i){
			continue;
		}
		reg = pleat->neighbors[i].region;
		if(reg.tag == GT_TwistTag){
			break;
		}
	}
	GT_TwistRegion *currTwist = reg.twist;
	if(GT_Polygon_contains(currTwist->n_sides, currTwist->vertices, a, 0)){
		return reg;
	}
	currTwist->scratch_data = !scratch;
	size_t cwi, ccwi;
	if(GT_Point_dot(GT_Point_sub(a, currTwist->vertices[0]), GT_Point_ccw(currTwist->neighbors[0].direction)) >= 0){
		//vertices[0] + t*directions[0] is a cw boundary of a around currTwist
		cwi = 0;
		for(size_t i = 1; i < currTwist->n_sides; ++i){
			if(GT_Point_dot(GT_Point_sub(a, currTwist->vertices[i]), GT_Point_ccw(currTwist->neighbors[i].direction)) >= 0){
				cwi = i;
			}else{
				ccwi = i;
				break;
			}
		}
	}else{
		//vertices[1] + t*directions[0] is a ccw boundary of a around currTwist
		ccwi = 0;
		for(size_t i = currTwist->n_sides - 1; i; --i){
			if(GT_Point_dot(GT_Point_sub(a, currTwist->vertices[i]), GT_Point_ccw(currTwist->neighbors[i].direction)) >= 0){
				cwi = i;
				break;
			}else{
				ccwi = i;
			}
		}
	}
	//so now finally we are searching at currTwist between cwi and ccwi
	//we can go along the pleat at cwi and search from the first twist node after [cwi].back_i
	//this will do a ccw spiral around a smaller and smaller set of nodes until the target is found
	GT_TaggedRegionPtr ret = findPaperRegion_rec(currTwist, cwi, ccwi, scratch, a);//TODO: if we build an explicit list of visited nodes we can tail call
	currTwist->scratch_data = scratch;
	external->scratch_data = scratch;
	return ret;
}

GT_TaggedRegionPtr findPaperRegion_rec(GT_TwistRegion *currTwist, size_t cwi, size_t ccwi, int scratch, GT_Point a){
	GT_TaggedRegionPtr region = currTwist->neighbors[cwi].region;
	GT_PleatRegion *pleat = region.pleat;
	if(GT_Polygon_contains(pleat->n_sides, pleat->vertices, a, 0)){
		return region;
	}
	size_t back_i = region.back_i, i;
	for(i = (back_i + 1)%pleat->n_sides; i != back_i; i = (i + 1)%pleat->n_sides){
		region = pleat->neighbors[i].region;
		if(region.tag == GT_TwistTag){
			break;
		}
	}
	GT_TwistRegion *nextTwist = region.twist;
	if(nextTwist->scratch_data != scratch){//loop detected
		return findPaperRegion_base(pleat, i, a);
	}
	size_t ccwi2 = region.back_i, cwi2 = ccwi2;
	if(GT_Polygon_contains(nextTwist->n_sides, nextTwist->vertices, a, 0)){
		return region;
	}
	nextTwist->scratch_data = !scratch;
	while(GT_Point_dot(GT_Point_sub(a, nextTwist->vertices[cwi2]), GT_Point_ccw(nextTwist->neighbors[cwi2].direction)) < 0){
		cwi2 = cwi2 ? cwi2 - 1 : nextTwist->n_sides - 1;
	}
	GT_TaggedRegionPtr ret = findPaperRegion_rec(nextTwist, cwi2, ccwi2, scratch, a);
	nextTwist->scratch_data = scratch;
	return ret;
}

GT_TaggedRegionPtr findPaperRegion_base(GT_PleatRegion *pleat, size_t i, GT_Point a){
	return pleat->neighbors[(i + 1)%pleat->n_sides].region;
}

GT_TwistRegion *GT_FlatRegion_spliceInTwist(GT_FlatRegion *flat_node, GT_TwistRegion *twist_node){
	//twist_node needs to be fully contained in flat_node (I might want a special function for this
	GT_Point *c_points = malloc((flat_node->n_sides + twist_node->n_sides)*sizeof(GT_Point));//TODO: check failed alloc
	size_t boundary_intersections = 0;
	GT_Polygon_intersectConvex(c_points, &boundary_intersections, flat_node->n_sides, flat_node->vertices, twist_node->n_sides, twist_node->vertices);
	free(c_points);
	GT_Point twist_center = GT_Polygon_interiorPoint(twist_node->n_sides, twist_node->vertices);
	//if the boundaries intersect, neither polygon can be contained in the other
	//if the boundaries do not intersect, if flat_node contains a point in the interior of twist_node then it contains all of twist_node
	//note that I believe the polygon intersection counts exactly the boundary crossings as boundary intersections which is what I want, but this may be wrong
	if(boundary_intersections || !GT_Polygon_contains(flat_node->n_sides, flat_node->vertices, twist_center, 0)){
		return NULL;
	}
	//now we have to split flat_node into twist_node->n_sides parts and add twist_node->n_sides pleat nodes as well as the actual twist node
	//but at each intersection of the new pleats with the existing borders of flat_node, what should be done?
	//for now I can just put in a triangular twist region and make a registry of these so they can easily be replaced with the proper twist that should go there
	//but if the pleat region hits more than one boundary pleat region, I will generate multiple triangular twist regions
	//if these overlap with existing synthetic triangular twists, they will be merged into synthetic rectangular twists
}

