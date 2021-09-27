![](gtess.png)
# GTess
GTess is a python library and pygame application for designing origami tesselations.
Currently, crease patterns must be generated manually using the API, as the application is
only capable of viewing them.

## Goal
The main goal of GTess is to be able to design a tesselated flower tower.
This is very frustrating to do by hand.
GTess will be able to help the user line up mutiple "units" within a crease pattern,
ensuring that their occupied paper regions do not overlap and their outgoing pleats dovetail.
It will then be able to report measurements of lines and angles and coordinates of points
to aid in actually folding the generated design.

Other features, such as creating "units" (eg the n sided flower tower as an intersection
of n split pleats) from simpler units or from scratch, specifying (or even automatically
generating) folded forms of units for a preview of the final tesselation, and so on,
may be added after the program is functional.



## Implementation Notes
The python code has extensive docstrings, but some general design notes follow here.

### Points and Transformations

GTess operates entirely on 2d data, because crease patterns are 2d by nature.
Automatically generating folded versions of crease patterns in general requires working in
3d, but crease pattern foldability in general is an NP complete problem so if it is
eventually implemented, some compromises would be required.

Because all data is 2d, points are stored as `x`, `y` coordinates and transformations are
stored as 2 by 3 matrices.  This is because GTess uses a common technique called
homogeneous coordinates to allow translating points to be accomplished via matrix
multiplication.  Instead of representing a point as `x`, `y`, we represent it as
`x`, `y`, `1`.  The redundant `1` coordinate is not actually stored, but having a
coordinate that is defined to be `1` is what allows us to encode translations as matrix
multiplications.  The bottom row of the matrix is always `[0, 0, 1]` so we omit it.

Be careful when editing the code: numpy arrays are gotcha factories.  The main gotchas
stem from the fact that getting a row or column of a matrix does not return a 2d array,
but rather a 1d array, which acts like a row vector, except that transposing it does not
make it a column vector.  `.reshape` is needed to make the 1d array row vector into a
column vector.

### Crease Pattern Representation

By definition a crease pattern is a planar graph with some additional constraints.
Each edge is etiher a raw edge or a crease with exterior dihedral angle between `0` and
`2\pi`.  For a flat-foldable crease pattern though, each crease has an exterior dihedral
angle of either `0` or `2\pi` for a mountain or valley fold respectively (assuming the
"exterior" of the model is on the back side compared to the crease pattern).

Therefore, we can store a crease pattern as a planar graph where each edge is
labeled as "mountain", "valley", or "raw".

Because we have a planar graph, embedded in euclidean space no less, we definitely want
to know about its regions.  The edges separate the plane into several regions.
Any region not bordering a raw edge must be convex because of additional constraints on
crease patterns.  However, we still may have concave regions depending on the shape of the
boundary, and in particular the boundless "exterior" region outside the raw edges
is always convex.  These graphs (or any planar 3-connected graph) are called polyhedral
graphs because if they are embedded in a sphere instead of a plane they represent the
edges, vertices, and faces of a polyhedron.

Each region has a list of constituent vertices in counterclockwise order around the
center of ther region, plust a list of edges in the same order, starting with the edge
from the first vertex to the second.

Each edge has a reference to the region on its left and right, as well as the index in
each region's edge list it occurs at (which is also the index in each's vertex list its
first terminus in counterclockwise order occurs at).  Edges are additionally labeled as
mountain, valley, or raw.

Vertices contain no extra information besides their position.  They are simply stored
in each region.  This results in multiple redundant copies of each vertex, but greatly
simplifies modifying regions.

### Determining Regions from Vertices and Edges

Obviously, the regions of a planar graph are determined by its vertices and edges, but
we want to store the regions explicitly since recomputing them all the time would be
very inefficient and make writing algorithms much harder.

There are two main algorithms for finding the regions in a planar graph (at least that
I came up with): the scanline algorithm and the breadth first search algorithm.

A scanline algorithm works by moving a conceptual vertical line across an object.
As we sweep over the object, in this case a graph embedded in the plane, we keep
track of some information, in this case regions, and update our information at each
vertex.  When a vertex has multiple outgoing edges on or ahead of the scanline,
new regions are created.  The top and bottom edges of each region are tracked.
When they eventually meet up again, the region is finalized.

However, there is a small problem with this algorithm.  Namely, if our graph has
concave regions, we may create two regions and later discover we have to merge them.
This is perfectly doable, but it would be better if we could avoid this altogether
with a different algorithm.

The BFS algorithm avoids needing to merge regions ever.  It first loops over the
list of undirected edges to generate a list of adjacent vertices for each vertex.
Then it sorts the neighbors of each vertex in counterclocwise order.
Regions can then be found simply by keeping track of what directed edges we have
visited.  We pick a starting vertex arbitrarily and then move along one of the
edges coming out of it arbitrarily.  Then we repeatedly go to the first point
in clockwise order around the current point from the previous point, until
we close the region.  This causes us to go around the region in the proper
counterclockwise order.  Then we pick the first unvisitied edge in clocwise order
around the current point, moving to another point with unvisited outgoing edges
if we have visited all of them.  Once we have visited all directed edges, we are
done.

Both of these algorithms are linear time and constant space, if we assume the the
number of edges of a region is bounded by a constant for the scanline algorithm
and the number of edges meeting at a vertex is bounded by a constant for the
BFS algorithm.  These assumptions are both mostly true for crease patterns.
Regions with many sides are more common than vertices with many edges though.

