#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "point.h"
#include "region.h"

int GT_checkTwistPleatAdj(const GT_TwistRegion *twist_region, size_t twist_i, const GT_PleatRegion *pleat_region, size_t pleat_i){
    if(twist_i >= twist_region->n_sides || pleat_i >= pleat_region->n_sides){
        return 0;
    }
    GT_TaggedRegionPtr region = twist_region->neighbors[twist_i].region;
    if(region.tag != GT_PleatTag || region.pleat != pleat_region || region.back_i != pleat_i){
        return 0;
    }
    region = pleat_region->neighbors[pleat_i].region;
    if(region.tag != GT_TwistTag || region.twist != twist_region || region.back_i != twist_i){
        return 0;
    }
    const GT_Point twist_a = twist_region->vertices[twist_i], twist_b = twist_region->vertices[(twist_i + 1)%twist_region->n_sides];
    if(GT_Point_sqdist(twist_a, pleat_region->vertices[(pleat_i + 1)%pleat_region->n_sides]) > GT_EPSILON
       || GT_Point_sqdist(twist_b, pleat_region->vertices[pleat_i]) > GT_EPSILON){
        return 0;
    }else if(fabs(GT_Point_cross(pleat_region->direction, twist_region->neighbors[twist_i].direction)) > GT_EPSILON){
        return 0;
    }
    return GT_checkProfilesMatch(twist_region->neighbors[twist_i].profile, pleat_region->neighbors[pleat_i].profile);
}

int GT_checkTwistFlatAdj(const GT_TwistRegion *twist_region, size_t twist_i, const GT_FlatRegion *flat_region, size_t flat_i){
    if(!twist_region->is_external || twist_i >= twist_region->n_sides || flat_i >= flat_region->n_sides){
        return 0;
    }
    GT_TaggedRegionPtr region = twist_region->neighbors[twist_i].region;
    if(region.tag != GT_FlatTag || region.flat != flat_region || region.back_i != flat_i){
        return 0;
    }
    region = flat_region->neighbors[flat_i].region;
    if(region.tag != GT_TwistTag || region.twist != twist_region || region.back_i != twist_i){
        return 0;
    }
    const GT_Point twist_a = twist_region->vertices[twist_i], twist_b = twist_region->vertices[(twist_i + 1)%twist_region->n_sides];
    if(GT_Point_sqdist(twist_a, flat_region->vertices[(flat_i + 1)%flat_region->n_sides]) > GT_EPSILON
       || GT_Point_sqdist(twist_b, flat_region->vertices[flat_i]) > GT_EPSILON){
        return 0;
    }
    return !twist_region->neighbors[twist_i].profile.n_pleats;
}

int GT_checkFlatPleatAdj(const GT_FlatRegion *flat_region, size_t flat_i, const GT_PleatRegion *pleat_region, size_t pleat_i){
    if(flat_i >= flat_region->n_sides || pleat_i >= pleat_region->n_sides){
        return 0;
    }
    GT_TaggedRegionPtr region = flat_region->neighbors[flat_i].region;
    if(region.tag != GT_PleatTag || region.pleat != pleat_region || region.back_i != pleat_i){
        return 0;
    }
    region = pleat_region->neighbors[pleat_i].region;
    if(region.tag != GT_FlatTag || region.flat != flat_region || region.back_i != flat_i){
        return 0;
    }
    const GT_Point flat_a = flat_region->vertices[flat_i], flat_b = flat_region->vertices[(flat_i + 1)%flat_region->n_sides];
    if(GT_Point_sqdist(flat_a, pleat_region->vertices[(pleat_i + 1)%pleat_region->n_sides]) > GT_EPSILON
       || GT_Point_sqdist(flat_b, pleat_region->vertices[pleat_i]) > GT_EPSILON){
        return 0;
    }
    return !pleat_region->neighbors[pleat_i].profile.n_pleats;
}

static inline int checkProfiles_advancePP0(int *on_pleat, size_t *i, size_t *j, const GT_PleatProfile **pleat_profile,
                                           double *coord, double *base, double *length,
                                           const GT_PleatRegion *self, GT_Point perp){
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
        if(++*i == self->n_sides){
            return 0;
        }else{
            *coord = *base = GT_Point_dot(self->vertices[*i], perp);
            *length = GT_Point_dot(self->vertices[(*i + 1)%self->n_sides], perp) - *coord;
            *pleat_profile = &self->neighbors[*i].profile;
            *on_pleat = 0;
        }
    }
    return 1;
}

static inline int checkProfiles_advancePP1(int *on_pleat, size_t *i, size_t *j, const GT_PleatProfile **pleat_profile,
                                           double *coord, double *base, double *length,
                                           const GT_PleatRegion *self, GT_Point perp){
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
            *coord = *base = GT_Point_dot(self->vertices[*i], perp);
            *length = GT_Point_dot(self->vertices[*i - 1], perp) - *coord;
            *pleat_profile = &self->neighbors[*i - 1].profile;
            *on_pleat = 0;
        }
    }
    return 1;
}

int GT_PleatRegion_checkProfiles(const GT_PleatRegion *self){
    if(fabs(GT_Point_cross(self->direction, GT_Point_sub(self->vertices[1], self->vertices[0]))) > GT_EPSILON){
        return 0;
    }
    if(self->neighbors[0].profile.n_pleats){
        return 0;
    }
    size_t p2i = 0;
    for(size_t i = 2; i + 1 < self->n_sides; ++i){
        if(fabs(GT_Point_cross(self->direction, GT_Point_sub(self->vertices[i + 1], self->vertices[i]))) < GT_EPSILON){
            p2i = i;
            break;
        }
    }
    if(!p2i || self->neighbors[p2i].profile.n_pleats){
        return 0;
    }
    /* We now have to iterate over the pleat profiles of the edges on the vertices[1] end of the PleatRegion
     * and also the vertices[0] end.  Within each pleat profile there are potentially many pleats, and we need to
     * check that every pleat lines up with either a pleat on the other end or a vertex on the other end.
     */
    GT_Point perp = GT_Point_cw(self->direction);
    double pp1_base = GT_Point_dot(self->vertices[p2i], perp);
    double pp0_base = GT_Point_dot(self->vertices[p2i + 1], perp);
    double pp1_coord = pp1_base;
    double pp0_coord = pp0_base;
    double pp1_length = GT_Point_dot(self->vertices[p2i - 1], perp) - pp1_coord;
    double pp0_length = GT_Point_dot(self->vertices[(p2i + 2)%self->n_sides], perp) - pp0_coord;
    size_t pp1_i = p2i, pp1_j = 0;
    size_t pp0_i = p2i + 1, pp0_j = 0;
    const GT_PleatProfile *pp1 = &self->neighbors[pp1_i - 1].profile, *pp0 = &self->neighbors[pp0_i].profile;
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
                                                self, perp);
        }
        if(advance_pp1){
            working &= checkProfiles_advancePP1(&pp1_on_pleat, &pp1_i, &pp1_j, &pp1, &pp1_coord, &pp1_base, &pp1_length,
                                                self, perp);
        }
    }
    if(pp1_i > 1){
        for(; pp1_i > 1; --pp1_i){
            if(self->neighbors[pp1_i - 1].profile.n_pleats){
                return 0;
            }
        }
    }else for(; pp0_i < self->n_sides; ++pp0_i){
        if(self->neighbors[pp0_i].profile.n_pleats){
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

static void GT_forEach_TwistRegion(GT_TwistRegion*, size_t back_i, int updated, void*, void (*)(GT_TaggedRegionPtr, void*));
static void GT_forEach_PleatRegion(GT_PleatRegion*, size_t back_i, int updated, void*, void (*)(GT_TaggedRegionPtr, void*));
static void GT_forEach_FlatRegion(GT_FlatRegion*, size_t back_i, int updated, void*, void (*)(GT_TaggedRegionPtr, void*));

static void GT_forEach_TwistRegion(GT_TwistRegion *node, size_t back_i, int updated, void *data, void (*callback)(GT_TaggedRegionPtr, void*)){
    if(node->scratch_data == updated){
        return;
    }
    node->scratch_data = updated;
    for(size_t i = 0; i < node->n_sides; ++i){
        GT_TaggedRegionPtr region = node->neighbors[i].region;
        GT_forEach_PleatRegion(region.pleat, region.back_i, updated, data, callback);
    }
    callback((GT_TaggedRegionPtr){.tag=GT_TwistTag, .twist=node, .back_i=back_i}, data);
}

static void GT_forEach_PleatRegion(GT_PleatRegion *node, size_t back_i, int updated, void *data, void (*callback)(GT_TaggedRegionPtr, void*)){
    if(node->scratch_data == updated){
        return;
    }
    node->scratch_data = updated;
    for(size_t i = 0; i < node->n_sides; ++i){
        GT_TaggedRegionPtr region = node->neighbors[i].region;
        if(region.tag == GT_TwistTag){
            GT_forEach_TwistRegion(region.twist, region.back_i, updated, data, callback);
        }else{
            GT_forEach_FlatRegion(region.flat, region.back_i, updated, data, callback);
        }
    }
    callback((GT_TaggedRegionPtr){.tag=GT_PleatTag, .pleat=node, .back_i=back_i}, data);
}

static void GT_forEach_FlatRegion(GT_FlatRegion *node, size_t back_i, int updated, void *data, void (*callback)(GT_TaggedRegionPtr, void*)){
    if(node->scratch_data == updated){
        return;
    }
    node->scratch_data = updated;
    for(size_t i = 0; i < node->n_sides; ++i){
        GT_TaggedRegionPtr region = node->neighbors[i].region;
        if(region.tag == GT_PleatTag){
			GT_forEach_PleatRegion(region.pleat, region.back_i, updated, data, callback);
		}//otherwise tnp refers to the exterior twist node which is already visited
    }
    callback((GT_TaggedRegionPtr){.tag=GT_FlatTag, .flat=node, .back_i=back_i}, data);
}

void GT_forEachPaperRegion(GT_TwistRegion *exterior, void *data, void (*callback)(GT_TaggedRegionPtr, void*)){
    int updated = exterior->scratch_data = !exterior->scratch_data;
    for(size_t i = 0; i < exterior->n_sides; ++i){
        GT_TaggedRegionPtr region = exterior->neighbors[i].region;
        if(region.tag == GT_PleatTag){
            GT_forEach_PleatRegion(region.pleat, region.back_i, updated, data, callback);
        }else{
            GT_forEach_FlatRegion(region.flat, region.back_i, updated, data, callback);
        }
    }
    callback((GT_TaggedRegionPtr){.tag=GT_TwistTag, .twist=exterior, .back_i=exterior->n_sides}, data);//meaningless back_i
}

static void deletePaper_addToCleanupList(GT_TaggedRegionPtr region, GT_TaggedRegionPtr *data){
    switch(region.tag){
        case GT_TwistTag:
			region.twist->neighbors[0].region = *data;
			*data = region;
            break;
        case GT_PleatTag:
			region.pleat->neighbors[0].region = *data;
			*data = region;
            break;
        case GT_FlatTag:
			region.flat->neighbors[0].region = *data;
			*data = region;
            break;
    }
}

void GT_deletePaper(GT_TwistRegion *exterior){
	GT_FlatRegion dummy = {.scratch_data= !exterior->scratch_data};
	GT_TaggedRegionPtr cleanup_head = {.tag=GT_FlatTag, .flat= &dummy};
    GT_forEachPaperRegion(exterior, &cleanup_head, (void(*)(GT_TaggedRegionPtr, void*))deletePaper_addToCleanupList);
    for(GT_TaggedRegionPtr next; cleanup_head.tag != GT_FlatTag || cleanup_head.flat != &dummy; cleanup_head = next){
		switch(cleanup_head.tag){
			case GT_TwistTag:
				next = cleanup_head.twist->neighbors[0].region;
				free(cleanup_head.twist);
				break;
			case GT_PleatTag:
				next = cleanup_head.pleat->neighbors[0].region;
				free(cleanup_head.pleat);
				break;
			case GT_FlatTag:
				next = cleanup_head.flat->neighbors[0].region;
				free(cleanup_head.flat);
				break;
		}
	}
}

