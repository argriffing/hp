// Functions downstream of the 2d hp bond graph.
//
// The main function in this module computes the densest k-subgraph
// of a 2d hp bond graph.
// This can be done efficiently because of the particular
// constraints on the graph (all of its vertices have degree less than three,
// except for at most two vertices which may be of degree three).
//
// Although this function is not of general interest,
// it depends on functions of more general interest,
// such as segmentation into connected components,
// calculation of graph girth, and the solutions of subset sum problems.
//
// The bsg_ function name prefix stands for "bond subgraph".

#include "connectedcomponents.h"
#include "subsetsum.h"
#include "graphgirth.h"

// Track properties of the connected components of the bond graph.
// These do not depend on vertex order.
typedef struct tagBSG_COMPONENT {
  int nvertices; // number of vertices in the connected component
  int nedges;    // number of edges in the connected component
  int ell;       // nedges - nvertices
  int girth;     // length of smallest cycle or -1 if no cycle exists
  int index;     // index into graph partition structure
} BSG_COMPONENT;


// The input argument k defines the number of vertices to be selected.
// The first k entries of the selection output array
// give the indices of the selected vertices.
// The number of edges in the induced subgraph will be written
// to the subgraph_edge_count argument.
int bsg_get_densest_subgraph(
    int k,
    int *selection_out, int *subgraph_edge_count
    )
{
}

