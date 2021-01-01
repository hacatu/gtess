#ifndef GTESS_PAPER_TRANSFORMERS_H
#define GTESS_PAPER_TRANSFORMERS_H
#include <stddef.h>
#include <stdio.h>
#include "region.h"

GT_TaggedRegionPtr GT_findPaperRegion(GT_TwistRegion *external, GT_Point a);
GT_TwistRegion *GT_FlatNode_spliceInTwist(GT_FlatRegion *flat_region, GT_TwistRegion *twist_region);

#endif //GTESS_PAPER_TRANSFORMERS_H

