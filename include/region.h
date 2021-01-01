#ifndef GTESS_REGION_H
#define GTESS_REGION_H
#include <stddef.h>
#include <stdio.h>
#include "point.h"

typedef struct{
    double mountain_fraction, width_fraction;
} GT_Pleat;

typedef struct{
    size_t n_pleats;
    GT_Pleat *pleats;
} GT_PleatProfile;

typedef struct GT_TwistNeighbor GT_TwistNeighbor;

typedef struct{
    size_t n_sides;
    int is_external:1, is_synthetic:1, scratch_data:1;
    GT_Point *vertices;
    GT_TwistNeighbor *neighbors;
} GT_TwistRegion;

/*
typedef struct{
    size_t twists_len;
    GT_TwistRegion *twists;
    size_t pleat_buf_len;
    GT_Pleat *pleat_buf;
} GT_TwistRegionBundle;
*/

typedef struct GT_PleatNeighbor GT_PleatNeighbor;

typedef struct{
    size_t n_sides;
    int scratch_data:1;
    GT_Point direction;
    GT_Point *vertices;
    GT_PleatNeighbor *neighbors;
} GT_PleatRegion;

typedef struct GT_FlatNeighbor GT_FlatNeighbor;

typedef struct{
    size_t n_sides;
    int scratch_data:1;
    GT_Point *vertices;
    GT_FlatNeighbor *neighbors;
} GT_FlatRegion;

typedef struct{
    enum{
        GT_TwistTag,
        GT_PleatTag,
        GT_FlatTag
    } tag;
    union{
        GT_TwistRegion *twist;
        GT_PleatRegion *pleat;
        GT_FlatRegion *flat;
    };
    size_t back_i;
} GT_TaggedRegionPtr;

struct GT_TwistNeighbor{
	GT_Point direction;
	GT_PleatProfile profile;
	GT_TaggedRegionPtr region;
};

struct GT_PleatNeighbor{
	GT_PleatProfile profile;
	GT_TaggedRegionPtr region;
};

struct GT_FlatNeighbor{
	GT_TaggedRegionPtr region;
};

///Execute a callback on every node in the paper, including the exterior node.  Current implementation is dfs with ccw node ordering.
void GT_forEachPaperRegion(GT_TwistRegion *exterior, void *data, void (*callback)(GT_TaggedRegionPtr region, void *data));
///Free all nodes and pleat arrays in the paper.  If nodes or pleat arrays were not allocated by malloc, this should not be used.
void GT_deletePaper(GT_TwistRegion *exterior);

///Check if a twist node and a pleat node are adjacent to each other at the specified indices and with matching pleat profiles
int GT_checkTwistPleatAdj(const GT_TwistRegion *twist_region, size_t twist_i, const GT_PleatRegion *pleat_region, size_t pleat_i);
///Check if a twist node and a flat node are adjacent to each other at the specified indices.  This is only allowed or possible if the twist node is the external node.
int GT_checkTwistFlatAdj(const GT_TwistRegion *twist_region, size_t twist_i, const GT_FlatRegion *flat_region, size_t flat_i);
///Check if a flat node and a pleat node are adjacent to each other at the specified indices and the pleat node does not have any pleats on the corresponding edge
int GT_checkFlatPleatAdj(const GT_FlatRegion *flat_region, size_t flat_i, const GT_PleatRegion *pleat_region, size_t pleat_i);
///Check if a pleat region's pleat profiles are consistent with a sequence of parallel pleats parallel to the direction of the pleat region
int GT_PleatRegion_checkProfiles(const GT_PleatRegion *self);
///Check if two pleat profiles match up (they will be corresponding not the same, ie one should be from a pleat region and the other from a twist region at coincident edges)
int GT_checkProfilesMatch(GT_PleatProfile a, GT_PleatProfile b);

#endif //GTESS_REGION_H

