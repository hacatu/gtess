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

typedef struct{
    size_t n_sides;
    int is_external:1, scratch_data:1;
    char data_buf[];
} GT_TwistRegion;

typedef struct{
    size_t twists_len;
    GT_TwistRegion *twists;
    size_t pleat_buf_len;
    GT_Pleat *pleat_buf;
} GT_TwistRegionBundle;

typedef struct{
    size_t n_sides;
    int scratch_data:1;
    GT_Point direction;
    char data_buf[];
} GT_PleatRegion;

typedef struct{
    size_t n_sides;
    int scratch_data:1;
    char data_buf[];
} GT_FlatRegion;

typedef GT_TwistRegion GT_TwistNode;
typedef GT_PleatRegion GT_PleatNode;
typedef GT_FlatRegion GT_FlatNode;

typedef struct{
    enum{
        GT_TwistTag,
        GT_PleatTag,
        GT_FlatTag
    } tag;
    union{
        GT_TwistNode *twist;
        GT_PleatNode *pleat;
        GT_FlatNode *flat;
    };
    size_t back_i;
} GT_TaggedNodePtr;

GT_Point *GT_TwistRegion_getVertices(GT_TwistRegion *self);
GT_Point *GT_TwistRegion_getDirections(GT_TwistRegion *self);
GT_PleatProfile *GT_TwistRegion_getPleatProfiles(GT_TwistRegion *self);
GT_TwistRegion *GT_TwistRegion_getNext(GT_TwistRegion *self);
GT_TaggedNodePtr *GT_TwistNode_getAdjPleats(GT_TwistNode *self);
GT_TwistNode *GT_TwistNode_getNext(GT_TwistNode *self);

GT_Point *GT_PleatRegion_getVertices(GT_PleatRegion *self);
GT_PleatProfile *GT_PleatRegion_getPleatProfiles(GT_PleatRegion *self);
GT_TaggedNodePtr *GT_PleatNode_getAdjNodes(GT_PleatNode *self);

GT_Point *GT_FlatRegion_getVertices(GT_FlatRegion *self);
GT_TaggedNodePtr *GT_FlatNode_getAdjPleats(GT_FlatNode *self);

void GT_forEachPaperNode(GT_TwistNode *exterior_node, void *data, void (*callback)(GT_TaggedNodePtr tnp, void *data));
void GT_deletePaper(GT_TwistNode *exterior_node);

int GT_checkTwistPleatAdj(const GT_TwistNode *twist_node, size_t twist_i, const GT_PleatNode *pleat_node, size_t pleat_i);
int GT_checkTwistFlatAdj(const GT_TwistNode *twist_node, size_t twist_i, const GT_FlatNode *flat_node, size_t flat_i);
int GT_checkFlatPleatAdj(const GT_FlatNode *flat_node, size_t flat_i, const GT_PleatNode *pleat_node, size_t pleat_i);
int GT_PleatRegion_checkProfiles(const GT_PleatRegion *self);
int GT_checkProfilesMatch(GT_PleatProfile a, GT_PleatProfile b);

GT_TwistNode *GT_makeExterior_10inSquare(void);

#endif //GTESS_REGION_H

