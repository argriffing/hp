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


#include "assert.h"

#include "connectedcomponents.h"
//#include "subsetsum.h"
#include "graphgirth.h"
#include "breadthfirst.h"

// Track properties of the connected components of the bond graph.
// These do not depend on vertex order.
/*
typedef struct tagBSG_COMPONENT {
  int nvertices; // number of vertices in the connected component
  int nedges;    // number of edges in the connected component
  int ell;       // nedges - nvertices
  int girth;     // length of smallest cycle or -1 if no cycle exists
  int index;     // index into graph partition structure
} BSG_COMPONENT;
*/

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
    int *va_trace, int *vb_trace,
    BFS_WS *bfs_ws, int *depth_ws
    )
{
  int nedges = ccgraph_get_component_nedges(ccgraph, component);
  int nvertices = ccgraph_get_component_nvertices(ccgraph, component);
  int ell = nedges - nvertices;
  int girth = -1;
  SUBGRAPH *subgraph = ccgraph->subgraph + component;

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
  int local_witness = -1;
  get_girth_and_vertex_conn(local_row_ptr, local_col_ind, nvertices,
      bfs_ws, depth_ws, &girth, &local_witness);
  int global_witness = subgraph->local_to_global[local_witness];

  // Put the entries of the smallest cycle into the local to global map.
  int ncycle = -1;
  get_smallest_cycle_ub(
      row_ptr, col_ind, global_witness,
      bfs_ws, depth_ws,
      va_trace, vb_trace,
      subgraph->local_to_global, &ncycle);
  assert(ncycle == girth);

  // Use depth first search to fill the rest of the connected component.
  int nfilled = bfs_fill(row_ptr, col_ind,
      bfs_ws->parent, subgraph->local_to_global, ncycle);
  assert(nfilled == nvertices);
  bfs_clear(bfs_ws->parent, subgraph->local_to_global, nfilled);

  int i;
  int v_local, v_global;
  int w_local, w_global;

  // Translate this global vertex ordering back into local coordinates.
  for (v_local=0; v_local<subgraph->nvertices; ++v_local) {
    v_global = subgraph->local_to_global[v_local];
    ccgraph->global_to_local[v_global] = v_local;
  }

  // Update the subgraph csr.
  int nedges_updated = local_row_ptr[0];
  for (v_local=0; v_local<subgraph->nvertices; ++v_local) {

    // Get the global index of the vertex.
    v_global = subgraph->local_to_global[v_local];

    // For each edge, add the local sink vertex to the compo col ind array.
    for (i=row_ptr[v_global]; i<row_ptr[v_global+1]; ++i) {
      w_global = col_ind[i];
      w_local = ccgraph->global_to_local[w_global];
      ccgraph->compo_col_ind[nedges_updated++] = w_local;
    }

    // Add to the compo row ptr list.
    local_row_ptr[v_local + 1] = nedges_updated;
  }

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
// Input workspaces:
// va_trace, vb_trace, bfs_ws, depth_ws
//
// Outputs:
// The first k entries of the selection output array selection_out
// give the selected vertices.
// Return the number of edges in the induced subgraph.
//
int bsg_get_densest_subgraph(
    const int *row_ptr, const int *col_ind, int nvertices,
    CCGRAPH *ccgraph,
    int k,
    BSG_COMPONENT *bsg_ws,
    int *va_trace, int *vb_trace,
    BFS_WS *bfs_ws, int *depth_ws,
    int *selection_out
    )
{
  // Initialize the selection to nonsensical values.
  int v;
  for (v=0; v<nvertices; ++v) {
    selection_out[v] = -1;
  }

  // If the number of requested vertices exceeds the total number of vertices
  // in the graph then return a negative number.
  if (k > nvertices) {
    return -1;
  }

  // Initialize the number of components and vertices
  // remaining for consideration.
  // Initialize the total score.
  int ncomponents_remaining = ccgraph->ncomponents;
  int nvertices_remaining = nvertices;
  int total_score = 0;

  // For each component, move the smallest cycle to the front of the array
  // and order all vertices so that each new vertex is adjacent
  // to at least one earlier vertex.
  // Also set some component attributes.
  int c;
  BSG_COMPONENT *bsg;
  for (c=0; c<ncomponents_remaining; ++c) {
    int computed_girth = _move_smallest_cycle_to_front(
        row_ptr, col_ind,
        ccgraph, c,
        va_trace, vb_trace,
        bfs_ws, depth_ws);
    bsg = bsg_ws + c;
    bsg->nvertices = ccgraph_get_component_nvertices(ccgraph, c);
    bsg->nedges = ccgraph_get_component_nedges(ccgraph, c);
    bsg->ell = bsg->nedges - bsg->nvertices;
    bsg->girth = computed_girth;
    bsg->index = c;
  }

  // Track the score.

  int v_local, v_global;
  int w_local, w_global;

  // Push isolated vertices to the back
  // while the number of available remaining vertices is less than k
  // and the index of the component under consideration is less
  // than the number of components remaining for consideration.
  component = 0;
  while (nvertices_remaining < k && component < ncomponents_remaining)
  {
    bsg = bsg_ws + component;
    if (bsg->nedges == 0) {
      BSG_COMPONENT tmp = *bsg;
      *bsg = bsg_ws[ncomponents_remaining-1];
      bsg_ws[ncomponents_remaining-1] = tmp;
      nvertices_remaining--;
      ncomponents_remaining--;
    } else {
      component++;
    }
  }

  // Look for the special component with more than one more edge than vertex.
  // If this component exists and contains at most k vertices
  // then we will add it to the selection.
  int special_component_index = -1;
  for (c=0; c<ncomponents; ++c) {
    bsg = bsg_ws + c;
    if (bsg->ell > 1) {
      assert(special_component_index == -1); // at most one is allowed
      special_component_index = c;
    }
  }
  if (special_component_index >= 0) {
    bsg = bsg_ws + special_component_index;
    // If the 
    if (bsg->nvertices <= k) {
      for (v_local=0; v_local<bsg->nvertices; ++v_local) {
        v_global = ccgraph_local_to_global(ccgraph, bsg->index, v_local);
      }
    }
  }

  // Move the special components towards the front.
  // There should be at most a single special component
  // with greater than one more edge than vertex.
  BSG_COMPONENT *special_compos;
  int nspecial_compos;
  for (component=0; component<ccgraph->ncomponents; ++component) {
    if (ell 


}

