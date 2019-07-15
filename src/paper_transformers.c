#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "paper_transformers.h"
#include "polygon.h"

static GT_TaggedNodePtr findPaperRegion_rec(GT_TwistNode *currTwist, size_t cwi, size_t ccwi, int scratch, GT_Point a);
static GT_TaggedNodePtr findPaperRegion_base(GT_PleatNode *pleat, size_t i, GT_Point a);

GT_TaggedNodePtr GT_findPaperRegion(GT_TwistNode *external_node, GT_Point a){
	GT_TaggedNodePtr *adjNodes = GT_TwistNode_getAdjPleats(external_node);
	GT_TaggedNodePtr tn;
	int scratch = external_node->scratch_data;
	external_node->scratch_data = !scratch;
	size_t i;
	for(i = 0; i < external_node->n_sides; ++i){
		tn = adjNodes[i];
		if(tn.tag == GT_PleatTag){
			break;
		}
	}
	if(i == external_node->n_sides){//paper is just one flat region
		return tn;
	}
	GT_PleatNode *pleat = tn.pleat;
	if(GT_Polygon_contains(pleat->n_sides, GT_PleatRegion_getVertices(pleat), a, 0)){
		return tn;
	}
	size_t external_i = tn.back_i;
	adjNodes = GT_PleatNode_getAdjNodes(pleat);
	for(size_t i = 0; i < pleat->n_sides; ++i){//technically we are guaranteed to find an internal twist node adjacent to pleat but check bounds anyway
		if(i == external_i){
			continue;
		}
		tn = adjNodes[i];
		if(tn.tag == GT_TwistTag){
			break;
		}
	}
	GT_TwistNode *currTwist = tn.twist;
	GT_Point *vertices = GT_TwistRegion_getVertices(currTwist);
	if(GT_Polygon_contains(currTwist->n_sides, vertices, a, 0)){
		return tn;
	}
	currTwist->scratch_data = !scratch;
	GT_Point *directions = GT_TwistRegion_getDirections(currTwist);
	size_t cwi, ccwi;
	if(GT_Point_dot(GT_Point_sub(a, vertices[0]), GT_Point_ccw(directions[0])) >= 0){
		//vertices[0] + t*directions[0] is a cw boundary of a around currTwist
		cwi = 0;
		for(size_t i = 1; i < currTwist->n_sides; ++i){
			if(GT_Point_dot(GT_Point_sub(a, vertices[i]), GT_Point_ccw(directions[i])) >= 0){
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
			if(GT_Point_dot(GT_Point_sub(a, vertices[i]), GT_Point_ccw(directions[i])) >= 0){
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
	GT_TaggedNodePtr ret = findPaperRegion_rec(currTwist, cwi, ccwi, scratch, a);//TODO: if we build an explicit list of visited nodes we can tail call
	currTwist->scratch_data = scratch;
	external_node->scratch_data = scratch;
	return ret;
}

GT_TaggedNodePtr findPaperRegion_rec(GT_TwistNode *currTwist, size_t cwi, size_t ccwi, int scratch, GT_Point a){
	GT_TaggedNodePtr tn = GT_TwistNode_getAdjPleats(currTwist)[cwi];
	GT_PleatNode *pleat = tn.pleat;
	if(GT_Polygon_contains(pleat->n_sides, GT_PleatRegion_getVertices(pleat), a, 0)){
		return tn;
	}
	size_t back_i = tn.back_i, i;
	for(i = (back_i + 1)%pleat->n_sides; i != back_i; i = (i + 1)%pleat->n_sides){
		tn = GT_PleatNode_getAdjNodes(pleat)[i];
		if(tn.tag == GT_TwistTag){
			break;
		}
	}
	GT_TwistNode *nextTwist = tn.twist;
	if(nextTwist->scratch_data != scratch){//loop detected
		return findPaperRegion_base(pleat, i, a);
	}
	size_t ccwi2 = tn.back_i, cwi2 = ccwi2;
	GT_Point *vertices = GT_TwistRegion_getVertices(nextTwist);
	if(GT_Polygon_contains(currTwist->n_sides, vertices, a, 0)){
		return tn;
	}
	nextTwist->scratch_data = !scratch;
	GT_Point *directions = GT_TwistRegion_getDirections(nextTwist);
	while(GT_Point_dot(GT_Point_sub(a, vertices[cwi2]), GT_Point_ccw(directions[cwi2])) < 0){
		cwi2 = cwi2 ? cwi2 - 1 : nextTwist->n_sides - 1;
	}
	GT_TaggedNodePtr ret = findPaperRegion_rec(nextTwist, cwi2, ccwi2, scratch, a);
	nextTwist->scratch_data = scratch;
	return ret;
}

GT_TaggedNodePtr findPaperRegion_base(GT_PleatNode *pleat, size_t i, GT_Point a){
	return GT_PleatNode_getAdjNodes(pleat)[(i + 1)%pleat->n_sides];
}

GT_TwistNode *GT_FlatNode_spliceInTwist(GT_FlatNode *flat_node, GT_TwistNode *twist_node);

