#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "point.h"
#include "region.h"

GT_Point *GT_TwistRegion_getVertices(GT_TwistRegion *self){
    return (GT_Point*)self->data_buf;
}

GT_Point *GT_TwistRegion_getDirections(GT_TwistRegion *self){
    return (GT_Point*)(self->data_buf + self->n_sides*sizeof(GT_Point));
}

GT_PleatProfile *GT_TwistRegion_getPleatProfiles(GT_TwistRegion *self){
    return (GT_PleatProfile*)(self->data_buf + 2*self->n_sides*sizeof(GT_Point));
}

GT_TwistRegion *GT_TwistRegion_getNext(GT_TwistRegion *self){
    return (GT_TwistRegion*)(self->data_buf + self->n_sides*(2*sizeof(GT_Point) + sizeof(GT_PleatProfile)));
}

GT_TaggedNodePtr *GT_TwistNode_getAdjPleats(GT_TwistNode *self){
    return (GT_TaggedNodePtr*)(self->data_buf + self->n_sides*(2*sizeof(GT_Point) + sizeof(GT_PleatProfile)));
}

GT_TwistNode *GT_TwistNode_getNext(GT_TwistNode *self){
    return (GT_TwistNode*)(self->data_buf + self->n_sides*(2*sizeof(GT_Point) + sizeof(GT_PleatProfile) + sizeof(GT_TaggedNodePtr)));
}

GT_Point *GT_PleatRegion_getVertices(GT_PleatRegion *self){
    return (GT_Point*)(self->data_buf);
}

GT_PleatProfile *GT_PleatRegion_getPleatProfiles(GT_PleatRegion *self){
    return (GT_PleatProfile*)(self->data_buf + self->n_sides*sizeof(GT_Point));
}

GT_TaggedNodePtr *GT_PleatNode_getAdjNodes(GT_PleatNode *self){
    return (GT_TaggedNodePtr*)(self->data_buf + self->n_sides*(sizeof(GT_Point) + sizeof(GT_PleatProfile)));
}

GT_Point *GT_FlatRegion_getVertices(GT_FlatRegion *self){
    return (GT_Point*)(self->data_buf);
}

GT_TaggedNodePtr *GT_FlatNode_getAdjPleats(GT_FlatNode *self){
    return (GT_TaggedNodePtr*)(self->data_buf + self->n_sides*sizeof(GT_Point));
}

int GT_checkTwistPleatAdj(const GT_TwistNode *twist_node, size_t twist_i, const GT_PleatNode *pleat_node, size_t pleat_i){
    if(twist_i >= twist_node->n_sides || pleat_i >= pleat_node->n_sides){
        return 0;
    }
    const GT_Point *twist_vertices = GT_TwistRegion_getVertices((GT_TwistRegion*)twist_node);
    const GT_Point *pleat_vertices = GT_PleatRegion_getVertices((GT_PleatRegion*)pleat_node);
    const GT_TaggedNodePtr *twist_pleats = GT_TwistNode_getAdjPleats((GT_TwistNode*)twist_node);
    const GT_TaggedNodePtr *pleat_adjs = GT_PleatNode_getAdjNodes((GT_PleatNode*)pleat_node);
    GT_TaggedNodePtr tnp = twist_pleats[twist_i];
    if(tnp.tag != GT_PleatTag || tnp.pleat != pleat_node || tnp.back_i != pleat_i){
        return 0;
    }
    tnp = pleat_adjs[pleat_i];
    if(tnp.tag != GT_TwistTag || tnp.twist != twist_node || tnp.back_i != twist_i){
        return 0;
    }
    const GT_Point twist_a = twist_vertices[twist_i], twist_b = twist_vertices[(twist_i + 1)%twist_node->n_sides];
    if(GT_Point_sqdist(twist_a, pleat_vertices[(pleat_i + 1)%pleat_node->n_sides]) > GT_EPSILON
       || GT_Point_sqdist(twist_b, pleat_vertices[pleat_i]) > GT_EPSILON){
        return 0;
    }else if(fabs(GT_Point_cross(pleat_node->direction, GT_TwistRegion_getDirections((GT_TwistRegion*)twist_node)[twist_i])) > GT_EPSILON){
        return 0;
    }
    return GT_checkProfilesMatch(GT_TwistRegion_getPleatProfiles((GT_TwistRegion*)twist_node)[twist_i],
                                 GT_PleatRegion_getPleatProfiles((GT_PleatRegion*)pleat_node)[pleat_i]);
}

int GT_checkTwistFlatAdj(const GT_TwistNode *twist_node, size_t twist_i, const GT_FlatNode *flat_node, size_t flat_i){
    if(!twist_node->is_external || twist_i >= twist_node->n_sides || flat_i >= flat_node->n_sides){
        return 0;
    }
    const GT_Point *twist_vertices = GT_TwistRegion_getVertices((GT_TwistRegion*)twist_node);
    const GT_Point *flat_vertices = GT_FlatRegion_getVertices((GT_FlatRegion*)flat_node);
    const GT_TaggedNodePtr *twist_pleats = GT_TwistNode_getAdjPleats((GT_TwistNode*)twist_node);
    const GT_TaggedNodePtr *flat_adjs = GT_FlatNode_getAdjPleats((GT_FlatNode*)flat_node);
    GT_TaggedNodePtr tnp = twist_pleats[twist_i];
    if(tnp.tag != GT_FlatTag || tnp.flat != flat_node || tnp.back_i != flat_i){
        return 0;
    }
    tnp = flat_adjs[flat_i];
    if(tnp.tag != GT_TwistTag || tnp.twist != twist_node || tnp.back_i != twist_i){
        return 0;
    }
    const GT_Point twist_a = twist_vertices[twist_i], twist_b = twist_vertices[(twist_i + 1)%twist_node->n_sides];
    if(GT_Point_sqdist(twist_a, flat_vertices[(flat_i + 1)%flat_node->n_sides]) > GT_EPSILON
       || GT_Point_sqdist(twist_b, flat_vertices[flat_i]) > GT_EPSILON){
        return 0;
    }
    return !GT_TwistRegion_getPleatProfiles((GT_TwistRegion*)twist_node)[twist_i].n_pleats;
}

int GT_checkFlatPleatAdj(const GT_FlatNode *flat_node, size_t flat_i, const GT_PleatNode *pleat_node, size_t pleat_i){
    if(flat_i >= flat_node->n_sides || pleat_i >= pleat_node->n_sides){
        return 0;
    }
    const GT_Point *flat_vertices = GT_FlatRegion_getVertices((GT_FlatRegion*)flat_node);
    const GT_Point *pleat_vertices = GT_PleatRegion_getVertices((GT_PleatRegion*)pleat_node);
    const GT_TaggedNodePtr *flat_pleats = GT_FlatNode_getAdjPleats((GT_FlatNode*)flat_node);
    const GT_TaggedNodePtr *pleat_adjs = GT_PleatNode_getAdjNodes((GT_PleatNode*)pleat_node);
    GT_TaggedNodePtr tnp = flat_pleats[flat_i];
    if(tnp.tag != GT_PleatTag || tnp.pleat != pleat_node || tnp.back_i != pleat_i){
        return 0;
    }
    tnp = pleat_adjs[pleat_i];
    if(tnp.tag != GT_FlatTag || tnp.flat != flat_node || tnp.back_i != flat_i){
        return 0;
    }
    const GT_Point flat_a = flat_vertices[flat_i], flat_b = flat_vertices[(flat_i + 1)%flat_node->n_sides];
    if(GT_Point_sqdist(flat_a, pleat_vertices[(pleat_i + 1)%pleat_node->n_sides]) > GT_EPSILON
       || GT_Point_sqdist(flat_b, pleat_vertices[pleat_i]) > GT_EPSILON){
        return 0;
    }
    return !GT_PleatRegion_getPleatProfiles((GT_PleatRegion*)pleat_node)[pleat_i].n_pleats;
}

static inline int checkProfiles_advancePP0(int *on_pleat, size_t *i, size_t *j, const GT_PleatProfile **pleat_profile,
                                           double *coord, double *base, double *length, size_t n,
                                           const GT_Point vertices[static n],
                                           const GT_PleatProfile pleat_profiles[static n], GT_Point perp){
    int advance = 1;
    if(*on_pleat){
        if(++*j < (*pleat_profile)->n_pleats){
            *coord = *base + (*pleat_profile)->pleats[*j].mountain_fraction*(*length);
            advance = 0;
        }
    }else if((*pleat_profile)->n_pleats){
        *j = 0;
        *coord = *base + (*pleat_profile)->pleats[*j].mountain_fraction*(*length);
        *on_pleat = 1;
        advance = 0;
    }
    if(advance){
        if(++*i == n){
            return 0;
        }else{
            *coord = *base = GT_Point_dot(vertices[*i], perp);
            *length = GT_Point_dot(vertices[(*i + 1)%n], perp) - *coord;
            *pleat_profile = pleat_profiles + *i;
            *on_pleat = 0;
        }
    }
    return 1;
}

static inline int checkProfiles_advancePP1(int *on_pleat, size_t *i, size_t *j, const GT_PleatProfile **pleat_profile,
                                           double *coord, double *base, double *length, size_t n,
                                           const GT_Point vertices[], const GT_PleatProfile pleat_profiles[],
                                           GT_Point perp){
    int advance = 1;
    if(*on_pleat){
        if((*j)--){
            *coord = *base + (*pleat_profile)->pleats[*j].mountain_fraction*(*length);
            advance = 0;
        }
    }else if((*pleat_profile)->n_pleats){
        *j = (*pleat_profile)->n_pleats - 1;
        *coord = *base + (*pleat_profile)->pleats[*j].mountain_fraction*(*length);
        *on_pleat = 1;
        advance = 0;
    }
    if(advance){
        if(--*i == 1){
            return 0;
        }else{
            *coord = *base = GT_Point_dot(vertices[*i], perp);
            *length = GT_Point_dot(vertices[*i - 1], perp) - *coord;
            *pleat_profile = pleat_profiles + *i - 1;
            *on_pleat = 0;
        }
    }
    return 1;
}

int GT_PleatRegion_checkProfiles(const GT_PleatRegion *self){
    const GT_Point *vertices = GT_PleatRegion_getVertices((GT_PleatRegion*)self);
    const GT_PleatProfile *pleat_profiles = GT_PleatRegion_getPleatProfiles((GT_PleatRegion*)self);
    if(fabs(GT_Point_cross(self->direction, GT_Point_sub(vertices[1], vertices[0]))) > GT_EPSILON){
        return 0;
    }
    if(pleat_profiles[0].n_pleats){
        return 0;
    }
    size_t p2i = 0;
    for(size_t i = 2; i + 1 < self->n_sides; ++i){
        if(fabs(GT_Point_cross(self->direction, GT_Point_sub(vertices[i + 1], vertices[i]))) < GT_EPSILON){
            p2i = i;
            break;
        }
    }
    if(!p2i || pleat_profiles[p2i].n_pleats){
        return 0;
    }
    /* We now have to iterate over the pleat profiles of the edges on the vertices[1] end of the PleatRegion
     * and also the vertices[0] end.  Within each pleat profile there are potentially many pleats, and we need to
     * check that every pleat lines up with either a pleat on the other end or a vertex on the other end.
     */
    GT_Point perp = GT_Point_cw(self->direction);
    double pp1_base = GT_Point_dot(vertices[p2i], perp);
    double pp0_base = GT_Point_dot(vertices[p2i + 1], perp);
    double pp1_coord = pp1_base;
    double pp0_coord = pp0_base;
    double pp1_length = GT_Point_dot(vertices[p2i - 1], perp) - pp1_coord;
    double pp0_length = GT_Point_dot(vertices[(p2i + 2)%self->n_sides], perp) - pp0_coord;
    size_t pp1_i = p2i, pp1_j = 0;
    size_t pp0_i = p2i + 1, pp0_j = 0;
    const GT_PleatProfile *pp1 = pleat_profiles + pp1_i - 1, *pp0 = pleat_profiles + pp0_i;
    int pp1_on_pleat = 0, pp0_on_pleat = 0;
    for(int working = 1; working;){
        int advance_pp1 = 0, advance_pp0 = 0;
        if(fabs(pp1_coord - pp0_coord) < GT_EPSILON){
            if(pp1_on_pleat && pp0_on_pleat){
                if(fabs(pp1->pleats[pp1_j].width_fraction + pp0->pleats[pp0_j].width_fraction) > GT_EPSILON){
                    return 0;
                }
            }
            advance_pp1 = advance_pp0 = 1;
        }else if(pp1_coord < pp0_coord){
            if(pp1_on_pleat){
                return 0;
            }
            advance_pp1 = 1;
        }else{
            if(pp0_on_pleat){
                return 0;
            }
            advance_pp0 = 1;
        }
        if(advance_pp0){
            working &= checkProfiles_advancePP0(&pp0_on_pleat, &pp0_i, &pp0_j, &pp0, &pp0_coord, &pp0_base, &pp0_length,
                                                self->n_sides, vertices, pleat_profiles, perp);
        }
        if(advance_pp1){
            working &= checkProfiles_advancePP1(&pp1_on_pleat, &pp1_i, &pp1_j, &pp1, &pp1_coord, &pp1_base, &pp1_length,
                                                self->n_sides, vertices, pleat_profiles, perp);
        }
    }
    if(pp1_i > 1){
        for(; pp1_i > 1; --pp1_i){
            if(pleat_profiles[pp1_i - 1].n_pleats){
                return 0;
            }
        }
    }else for(; pp0_i < self->n_sides; ++pp0_i){
        if(pleat_profiles[pp0_i].n_pleats){
            return 0;
        }
    }
    return 1;
}

int GT_checkProfilesMatch(GT_PleatProfile a, GT_PleatProfile b){
    if(a.n_pleats != b.n_pleats){
        return 0;
    }else for(size_t i = 0; i < a.n_pleats; ++i){
        if(fabs(a.pleats[i].mountain_fraction + b.pleats[a.n_pleats - 1 - i].mountain_fraction - 1) > GT_EPSILON){
            return 0;
        }else if(fabs(a.pleats[i].width_fraction + b.pleats[a.n_pleats - 1 - i].width_fraction) > GT_EPSILON){
            return 0;
        }
    }
    return 1;
}

static void GT_forEach_TwistNode(GT_TwistNode*, size_t back_i, int updated, void*, void (*)(GT_TaggedNodePtr, void*));
static void GT_forEach_PleatNode(GT_PleatNode*, size_t back_i, int updated, void*, void (*)(GT_TaggedNodePtr, void*));
static void GT_forEach_FlatNode(GT_FlatNode*, size_t back_i, int updated, void*, void (*)(GT_TaggedNodePtr, void*));

static void GT_forEach_TwistNode(GT_TwistNode *node, size_t back_i, int updated, void *data, void (*callback)(GT_TaggedNodePtr, void*)){
    if(node->scratch_data == updated){
        return;
    }
    node->scratch_data = updated;
    GT_TaggedNodePtr *adjs = GT_TwistNode_getAdjPleats(node);
    for(size_t i = 0; i < node->n_sides; ++i){
        GT_TaggedNodePtr tnp = adjs[i];
        GT_forEach_PleatNode(tnp.pleat, tnp.back_i, updated, data, callback);
    }
    callback((GT_TaggedNodePtr){.tag=GT_TwistTag, .twist=node, .back_i=back_i}, data);
}

static void GT_forEach_PleatNode(GT_PleatNode *node, size_t back_i, int updated, void *data, void (*callback)(GT_TaggedNodePtr, void*)){
    if(node->scratch_data == updated){
        return;
    }
    node->scratch_data = updated;
    GT_TaggedNodePtr *adjs = GT_PleatNode_getAdjNodes(node);
    for(size_t i = 0; i < node->n_sides; ++i){
        GT_TaggedNodePtr tnp = adjs[i];
        if(tnp.tag == GT_TwistTag){
            GT_forEach_TwistNode(tnp.twist, tnp.back_i, updated, data, callback);
        }else{
            GT_forEach_FlatNode(tnp.flat, tnp.back_i, updated, data, callback);
        }
    }
    callback((GT_TaggedNodePtr){.tag=GT_PleatTag, .pleat=node, .back_i=back_i}, data);
}

static void GT_forEach_FlatNode(GT_FlatNode *node, size_t back_i, int updated, void *data, void (*callback)(GT_TaggedNodePtr, void*)){
    if(node->scratch_data == updated){
        return;
    }
    node->scratch_data = updated;
    GT_TaggedNodePtr *adjs = GT_FlatNode_getAdjPleats(node);
    for(size_t i = 0; i < node->n_sides; ++i){
        GT_TaggedNodePtr tnp = adjs[i];
        if(tnp.tag == GT_PleatTag){
			GT_forEach_PleatNode(tnp.pleat, tnp.back_i, updated, data, callback);
		}//otherwise tnp refers to the exterior twist node which is already visited
    }
    callback((GT_TaggedNodePtr){.tag=GT_FlatTag, .flat=node, .back_i=back_i}, data);
}

void GT_forEachPaperNode(GT_TwistNode *exterior_node, void *data, void (*callback)(GT_TaggedNodePtr, void*)){
    GT_TaggedNodePtr *exterior_adjs = GT_TwistNode_getAdjPleats(exterior_node);
    int updated = exterior_node->scratch_data = !exterior_node->scratch_data;
    for(size_t i = 0; i < exterior_node->n_sides; ++i){
        GT_TaggedNodePtr tnp = exterior_adjs[i];
        if(tnp.tag == GT_PleatTag){
            GT_forEach_PleatNode(tnp.pleat, tnp.back_i, updated, data, callback);
        }else{
            GT_forEach_FlatNode(tnp.flat, tnp.back_i, updated, data, callback);
        }
    }
    callback((GT_TaggedNodePtr){.tag=GT_TwistTag, .twist=exterior_node, .back_i=exterior_node->n_sides}, data);
}

static void deletePaper_addToCleanupList(GT_TaggedNodePtr tnp, GT_TaggedNodePtr *data){
    switch(tnp.tag){
        case GT_TwistTag:
			GT_TwistNode_getAdjPleats(tnp.twist)[0] = *data;
			*data = tnp;
            break;
        case GT_PleatTag:
			GT_PleatNode_getAdjNodes(tnp.pleat)[0] = *data;
			*data = tnp;
            break;
        case GT_FlatTag:
			GT_FlatNode_getAdjPleats(tnp.flat)[0] = *data;
			*data = tnp;
            break;
    }
}

void GT_deletePaper(GT_TwistNode *exterior_node){
	GT_FlatNode dummy = {.scratch_data= !exterior_node->scratch_data};
	GT_TaggedNodePtr cleanup_head = {.tag=GT_FlatTag, .flat= &dummy};
    GT_forEachPaperNode(exterior_node, &cleanup_head, (void(*)(GT_TaggedNodePtr, void*))deletePaper_addToCleanupList);
    for(GT_TaggedNodePtr next; cleanup_head.tag != GT_FlatTag || cleanup_head.flat != &dummy; cleanup_head = next){
		switch(cleanup_head.tag){
			case GT_TwistTag:
				next = GT_TwistNode_getAdjPleats(cleanup_head.twist)[0];
				free(cleanup_head.twist);
				break;
			case GT_PleatTag:
				next = GT_PleatNode_getAdjNodes(cleanup_head.pleat)[0];
				free(cleanup_head.pleat);
				break;
			case GT_FlatTag:
				next = GT_FlatNode_getAdjPleats(cleanup_head.flat)[0];
				free(cleanup_head.flat);
				break;
		}
	}
}

