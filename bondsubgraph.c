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
#include "breadthfirst.h"

// Track properties of the connected components of the bond graph.
// These do not depend on vertex order.
typedef struct tagBSG_COMPONENT {
  int nvertices; // number of vertices in the connected component
  int nedges;    // number of edges in the connected component
  int ell;       // nedges - nvertices
  int girth;     // length of smallest cycle or -1 if no cycle exists
  int index;     // index into graph partition structure
} BSG_COMPONENT;


// Reorder vertices within a connected component of an undirected graph
// so that the vertices in the smallest cycle are before the other vertices,
// and so that each vertex is adjacent to a vertex somewhere
// before it in the vertex list.
// The length of the smallest cycle, also called the girth, is returned.
// If the component has no cycle, then return -1 to represent an infinite girth
// and do not change the vertex order.
int _move_smallest_cycle_to_front(
    const int *row_ptr, const int *col_ind,
    CCGRAPH *ccgraph, int component,
    BFS_WS *bfs_ws, int *depth_ws
    )
{
  int nedges = ccgraph_get_component_nedges(ccgraph, component);
  int nvertices = ccgraph_get_component_nvertices(ccgraph, component);
  int ell = nedges - nvertices;
  int girth = -1;

  // If the number of vertices is much greater than the number
  // of edges, then we do not have a connected component.
  assert(ell >= -1);

  // If the number of vertices is greater than the number of edges
  // then we have something like a point or path or tree.
  // In any case we have no cycles, so do not change the order,
  // and return -1 indicating infinite girth.
  if (ell == -1) {
    girth = -1;
    return girth;
  }

  // If the number of vertices is equal to the number of edges
  // then we have a pure cycle.
  // Because there are no other vertices in this component,
  // this cycle is already at the front of the vertex list,
  // so we do not have to reorder any vertices.
  // The girth is the number of vertices in this cycle.
  if (ell == 0) {
    girth = nvertices;
    return girth;
  }

  // For more interesting connected components,
  // each cycle includes at least one vertex of degree at least three.
  // So we can use a modified topo sort to check each of these vertices
  // to see which one is part of the smallest cycle in the component.
  int *local_row_ptr = ccgraph_get_component_row_ptr(ccgraph, component);
  int *local_col_ind = ccgraph_get_component_col_ind(ccgraph, component);
  get_girth_and_vertex_conn(local_row_ptr, local_col_ind, nvertices,
      bfs_ws, depth_ws, &girth, &local_witness);
  int global_witness = ccgraph->local_to_global[local_witness];
  get_smallest_cycle_ub(
      row_ptr, col_ind, int nvertices, int r,
      BFS_WS *bfs_ws, int *depth_ws,
      int *cycle_out, int *ncycle_out
      );

  
  nvertices = bfs_fill(const int *row_ptr, const int *col_ind,
      int *parent, int *history, int nseeds);

  return girth;
}


// Get the set of vertices defining the densest induced k-subgraph.
//
// Inputs:
// A graph in csr format.
// A decomposition of the graph into connected components.
// An argument k defining the number of vertices to be selected.
// Another input is a workspace for bsg components;
// this is just a BSG_COMPONENT array at least as long as the
// number of components (which is bounded above by the number of vertices).
//
// Outputs:
// The first k entries of the selection output array selection_out
// give the indices of the selected vertices.
// The number of edges in the induced subgraph will be written
// to the subgraph_edge_count argument.
//
int bsg_get_densest_subgraph(
    const int *row_ptr, const int *col_ind, int nvertices,
    CCGRAPH *pccgraph,
    int k,
    BSG_COMPONENT *bsg_ws, // bond graph component workspace
    int *selection_out, int *subgraph_edge_count
    )
{
  int i;
  for (i=0; i<pccgraph->ncomponents; ++i) {
    bsg_ws;
  }
}

