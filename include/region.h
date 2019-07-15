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
    /* GT_Point vertices[n_sides];
     * GT_Point directions[n_sides];
     * GT_PleatProfile pleatProfiles[n_sides];
     * //GT_TwistNode adds: GT_TaggedNodePtr adjPleats[n_sides];
     */
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
    /* GT_Point vertices[n_sides];
     * GT_PleatProfile pleatProfiles[n_sides];
     * //GT_PleatNode adds: GT_TaggedNodePtr adjNodes[n_sides];
     */
} GT_PleatRegion;

typedef struct{
    size_t n_sides;
    int scratch_data:1;
    char data_buf[];
    /* GT_Point vertices[n_sides];
     * //GT_FlatNode adds: GT_TaggedNodePtr adjPleats[n_sides];
     */
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

///Get a pointer to the beginning of the array of vertices, stored in the data buffer
GT_Point *GT_TwistRegion_getVertices(GT_TwistRegion *self);
///Get a pointer to the beginning of the array of directions, stored in the data buffer
GT_Point *GT_TwistRegion_getDirections(GT_TwistRegion *self);
///Get a pointer to the beginning of the array of pleat profiles, stored in the data buffer
GT_PleatProfile *GT_TwistRegion_getPleatProfiles(GT_TwistRegion *self);
///Compute the address of the next twist region structure in a buffer of twist regions since they have variable size
GT_TwistRegion *GT_TwistRegion_getNext(GT_TwistRegion *self);
///Get a pointer to the beginning of the array of tagged pointers to adjacent nodes (which are all pleat nodes unless self is external), stored in the data buffer (not present in GT_TwistRegions)
GT_TaggedNodePtr *GT_TwistNode_getAdjPleats(GT_TwistNode *self);
///Compute the address of the next twist node structure in a buffer of twist nodes since they have variable size.  Twist nodes are twist regions that also have adjacent node information
GT_TwistNode *GT_TwistNode_getNext(GT_TwistNode *self);

///Get a pointer to the beginning of the array of vertices, stored in the data buffer
GT_Point *GT_PleatRegion_getVertices(GT_PleatRegion *self);
/**
 * Get a pointer to the beginning of the array of pleat profiles, stored in the data buffer
 * 
 * There must be at least two sides which are parallel to the direction of the pleat and have no pleats on them
 */
GT_PleatProfile *GT_PleatRegion_getPleatProfiles(GT_PleatRegion *self);
/**
 * Get a pointer to the beginning of the array of tagged pointers to adjacent nodes, stored in the data buffer (not present in GT_PleatRegions)
 * 
 * The nodes adjacent to the two parallel sides with no pleats should be flat regions and the other adjacent nodes should be twist nodes
 */
GT_TaggedNodePtr *GT_PleatNode_getAdjNodes(GT_PleatNode *self);

///Get a pointer to the beginning of the array of vertices, stored in the data buffer
GT_Point *GT_FlatRegion_getVertices(GT_FlatRegion *self);
///Get a pointer to the beginning of the array of tagged pointers to adjacent nodes (which are all pleat nodes or the external node), stored in the data buffer (not present in GT_FlatRegions)
GT_TaggedNodePtr *GT_FlatNode_getAdjPleats(GT_FlatNode *self);

///Execute a callback on every node in the paper, including the exterior node.  Current implementation is dfs with ccw node ordering.
void GT_forEachPaperNode(GT_TwistNode *exterior_node, void *data, void (*callback)(GT_TaggedNodePtr tnp, void *data));
///Free all nodes and pleat arrays in the paper.  If nodes or pleat arrays were not allocated by malloc, this should not be used.
void GT_deletePaper(GT_TwistNode *exterior_node);

///Check if a twist node and a pleat node are adjacent to each other at the specified indices and with matching pleat profiles
int GT_checkTwistPleatAdj(const GT_TwistNode *twist_node, size_t twist_i, const GT_PleatNode *pleat_node, size_t pleat_i);
///Check if a twist node and a flat node are adjacent to each other at the specified indices.  This is only allowed or possible if the twist node is the external node.
int GT_checkTwistFlatAdj(const GT_TwistNode *twist_node, size_t twist_i, const GT_FlatNode *flat_node, size_t flat_i);
///Check if a flat node and a pleat node are adjacent to each other at the specified indices and the pleat node does not have any pleats on the corresponding edge
int GT_checkFlatPleatAdj(const GT_FlatNode *flat_node, size_t flat_i, const GT_PleatNode *pleat_node, size_t pleat_i);
///Check if a pleat region's pleat profiles are consistent with a sequence of parallel pleats parallel to the direction of the pleat region
int GT_PleatRegion_checkProfiles(const GT_PleatRegion *self);
///Check if two pleat profiles match up (they will be corresponding not the same, ie one should be from a pleat region and the other from a twist region at coincident edges)
int GT_checkProfilesMatch(GT_PleatProfile a, GT_PleatProfile b);

//GT_TwistNode *GT_makeExterior_10inSquare(void);

#endif //GTESS_REGION_H

