#include "stdlib.h"
#include "stdio.h"
#include "assert.h"

#include "connectedcomponents.h"
#include "breadthfirst.h"


void _ccgraph_check_component(CCGRAPH *p, int component) {
  assert(0 <= component);
  assert(component < p->ncomponents);
}

int ccgraph_get_component_nvertices(CCGRAPH *p, int component) {
  _ccgraph_check_component(p, component);
  return p->subgraph[component].nvertices;
}

int ccgraph_get_component_nedges(CCGRAPH *p, int component) {
  _ccgraph_check_component(p, component);
  int nvertices = p->subgraph[component].nvertices;
  int *row_ptr = p->subgraph[component].row_ptr;
  return row_ptr[nvertices] - row_ptr[0];
}

int *ccgraph_get_component_row_ptr(CCGRAPH *p, int component) {
  _ccgraph_check_component(p, component);
  return p->subgraph[component].row_ptr;
}

int *ccgraph_get_component_col_ind(CCGRAPH *p, int component) {
  _ccgraph_check_component(p, component);
  return p->compo_col_ind;
}

int ccgraph_local_to_global(CCGRAPH *p, int component, int v_local) {
  _ccgraph_check_component(p, component);
  return p->subgraph[component].local_to_global[v_local];
}


// This function needs to be called only once per application,
// if the max number of vertices and max number of edges is known.
// You do not need to call this function every time you want
// to compute the connected components.
void ccgraph_init(CCGRAPH *p, int max_nvertices, int max_nedges)
{
  p->subgraph = (SUBGRAPH *) malloc(max_nvertices * sizeof(CCGRAPH));
  p->global_to_local = (int *) malloc(max_nvertices * sizeof(int));
  p->local_to_global = (int *) malloc(max_nvertices * sizeof(int));
  p->compo_row_ptr = (int *) malloc((max_nvertices + 1) * sizeof(int));
  p->compo_col_ind = (int *) malloc(max_nedges * sizeof(int));
}

void ccgraph_destroy(CCGRAPH *p)
{
  free(p->subgraph);
  free(p->global_to_local);
  free(p->local_to_global);
  free(p->compo_row_ptr);
  free(p->compo_col_ind);
}


// This is a helper function for flood filling a single component.
// The function name underscore prefix indicates that this function is not
// part of the public interface of this module.
// C may have keywords to enforce this somehow but I'm not doing that.
//
// The parent vertices are assumed to have been set to -1,
// and all vertices in the connected component of the root are assumed to be -1
// representing membership in a currently unknown component.
// The set of vertices in the component is also tracked.
// The local to global map and the global to local map should each
// have -1 entries to represent unknown values.
//
void _ccgraph_compute_component(CCGRAPH *p,
    const int *row_ptr, const int *col_ind,
    BFS_WS *bfs_ws,
    int root, int *component_labels
    )
{
  int i;
  int v_local, v_global;
  int w_local, w_global;

  // The current label is the number of components in the ccgraph.
  int label = p->ncomponents;

  // Pick the subgraph that corresponds to this component.
  SUBGRAPH *subgraph = p->subgraph + p->ncomponents;

  // Set up the subgraph indirection relative to the ccgraph object.
  subgraph->local_to_global = p->local_to_global + p->nvertices;
  subgraph->row_ptr = p->compo_row_ptr + p->nvertices;

  // Compute the connected component.
  subgraph->local_to_global[0] = root;
  subgraph->nvertices = 1;
  subgraph->nvertices = bfs_fill(row_ptr, col_ind,
      bfs_ws->parent, subgraph->local_to_global, subgraph->nvertices);
  
  // Set the labels.
  for (v_local=0; v_local<subgraph->nvertices; ++v_local) {
    v_global = subgraph->local_to_global[v_local];
    assert(component_labels[v_global] == -1);
    component_labels[v_global] = label;
  }
  
  // Reset parent vertices.
  bfs_clear(bfs_ws->parent, subgraph->local_to_global, subgraph->nvertices);

  // Define the map from global index to local index.
  for (v_local=0; v_local<subgraph->nvertices; ++v_local) {
    v_global = subgraph->local_to_global[v_local];
    p->global_to_local[v_global] = v_local;
  }

  // Update the subgraph csr.
  for (v_local=0; v_local<subgraph->nvertices; ++v_local) {

    // Get the global index of the vertex.
    v_global = subgraph->local_to_global[v_local];

    // For each edge, add the local sink vertex to the compo col ind array.
    for (i=row_ptr[v_global]; i<row_ptr[v_global+1]; ++i) {
      w_global = col_ind[i];
      w_local = p->global_to_local[w_global];
      printf("number of edges: %d\n", p->nedges);
      p->compo_col_ind[p->nedges++] = w_local;
    }

    // Add to the compo row ptr list.
    p->compo_row_ptr[p->nvertices + 1] = p->nedges;

    // Increment the number of vertices for which the edges are tracked.
    p->nvertices++;
  }

  // Increment the number of components for which the local graph is available.
  p->ncomponents++;
}


// This function assumes that the ccgraph has already been initialized.
void ccgraph_compute(CCGRAPH *p,
    const int *row_ptr, const int *col_ind, int nvertices,
    BFS_WS *bfs_ws,
    int *component_labels
    )
{
  p->ncomponents = 0;
  p->nvertices = 0;
  p->nedges = 0;
  p->compo_row_ptr[p->nvertices] = 0;

  int v;

  // Set all component labels and all parent vertices to -1.
  // Set the local/global index conversion map entries to -1.
  for (v=0; v<nvertices; ++v) {
    component_labels[v] = -1;
    bfs_ws->parent[v] = -1;
    p->local_to_global[v] = -1;
    p->global_to_local[v] = -1;
  }

  // For each vertex that has not been assigned a component,
  // flood fill a component rooted at that vertex.
  // The helper function does post-processing including
  // incrementing the number of components.
  int root;
  for (root=0; root<nvertices; ++root) {
    if (component_labels[root] < 0) {
      _ccgraph_compute_component(p, row_ptr, col_ind,
          bfs_ws,
          root, component_labels);
    }
  }
}

