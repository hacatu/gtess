#ifndef GTESS_PAPER_TRANSFORMERS_H
#define GTESS_PAPER_TRANSFORMERS_H
#include <stddef.h>
#include <stdio.h>
#include "point.h"

GT_TaggedNodePtr GT_findPaperRegion(GT_Point *a);
GT_TwistNode *GT_FlatNode_spliceInTwist(GT_FlatNode *flat_node, GT_TwistNode *twist_node);

#endif //GTESS_PAPER_TRANSFORMERS_H

