#ifndef GTESS_QTREE_H
#define GTESS_QTREE_H

#include "iteractive.h"

typedef struct GT_Qtree GT_Qtree;
struct GT_Qtree{
	double mid_x, mid_y;
	GT_Qtree children[2][2];
	GT_Hitrect *hitrects;
	size_t len, cap, total_hitrect_n;
};

//TODO: quadtrees are useful, but currently they are not needed: I have an algorithm for finding
//what paper region was clicked when that is relevant and for popup menus and stuff there will
//be few enough buttons and things will be happening in user time so speed won't be a concern.
//Quadtrees will be properly implemented later if necessary.

#endif

