// Conventions:
// Component label -1 means unlabeled.
// Parent -2 means root.
// Parent -1 means not visited by depth first search.
//
// Connected components in csr format.

#ifndef CONNECTEDCOMPONENTS_HEADER
#define CONNECTEDCOMPONENTS_HEADER

#include "breadthfirst.h"


// The idea of this structure is to avoid malloc/free in an inner loop.
// This structure should not own the memory referenced by any of its pointers.
typedef struct tagSUBGRAPH {

  // This is the number of local vertices.
  int nvertices;

  // Points somewhere in ccgraph.compo_row_ptr.
  int *row_ptr;

  // Points somewhere in ccgraph.local_to_global.
  // Maps local vertex indices to global vertex indices.
  int *local_to_global;

} SUBGRAPH;


typedef struct tagCCGRAPH {
  int nvertices;
  int ncomponents;
  int nedges;
  SUBGRAPH *subgraph; // length is nvertices (as upper bound on ncomponents)
  int *global_to_local; // length is nvertices
  int *local_to_global; // length is nvertices
  int *compo_row_ptr; // length is nvertices + 1
  int *compo_col_ind; // length is nedges
} CCGRAPH;

int ccgraph_get_component_nvertices(CCGRAPH *p, int component);

void ccgraph_get_component_degree_min_max(CCGRAPH *p, int component,
    int *min_degree, int *max_degree);

int ccgraph_get_component_nedges(CCGRAPH *p, int component);

int ccgraph_get_component_nedges_undirected(CCGRAPH *p, int component);

int *ccgraph_get_component_row_ptr(CCGRAPH *p, int component);

int *ccgraph_get_component_col_ind(CCGRAPH *p, int component);

int ccgraph_local_to_global(CCGRAPH *p, int component, int v_local);


// This function needs to be called only once per application,
// if the max number of vertices and max number of edges is known.
// You do not need to call this function every time you want
// to compute the connected components.
void ccgraph_init(CCGRAPH *p, int max_nvertices, int max_nedges);

void ccgraph_destroy(CCGRAPH *p);

// This function assumes that the ccgraph has already been initialized.
void ccgraph_compute(CCGRAPH *p,
    const int *row_ptr, const int *col_ind, int nvertices,
    BFS_WS *bfs_ws,
    int *component_labels
    );


#endif

