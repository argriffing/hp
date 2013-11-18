#include "stdio.h"
#include "stdlib.h"

#include "sparsetools.h"
#include "connectedcomponents.h"

int t0()
{
  /*
   *   0        4--5--6
   *                           7
   *   1--2         3         / \
   *                         /   \
   *     8                  11----9
   *    / \
   *   /   \
   * 13-----12
   *   \   /
   *    \ /
   *     10
   *
   * */
  int lil[] = {
    1, 2,
    4, 5,
    5, 6,
    7, 9,
    9, 11,
    11, 7,
    8, 12,
    12, 10,
    10, 13,
    13, 8,
    12, 13,
  };
  int nvertices = 14;
  int nedges = 11;

  // Allocate the arrays that define the graph.
  int *row_ptr = (int *) malloc((nvertices + 1) * sizeof(int));
  int *col_ind = (int *) malloc(2 * nedges * sizeof(int));

  // Make a csr format graph.
  lil_to_csr(nvertices, nedges, lil, row_ptr, col_ind);

  // Allocate things for the connected components search.
  int *parent_ws = (int *) malloc(nvertices * sizeof(int));
  int *deck_ws = (int *) malloc(nvertices * sizeof(int));
  int *next_ws = (int *) malloc(nvertices * sizeof(int));
  int *component_labels = (int *) malloc(nvertices * sizeof(int));

  // Init the connected components search object.
  CCGRAPH g;
  ccgraph_init(&g, nvertices, 2 * nedges);

  // Do the connected components search.
  ccgraph_compute(&g, row_ptr, col_ind, nvertices,
      parent_ws, deck_ws, next_ws,
      component_labels);

  printf("%d connected components were detected\n", g.ncomponents);
  printf("\n");

  printf("connected component graphs:\n");
  int component, i;
  int v_local, v_global;
  int w_local, w_global;
  for (component=0; component<g.ncomponents; ++component) {
    printf("component %d:\n", component);
    int compo_nvertices = ccgraph_get_component_nvertices(&g, component);
    int *compo_row_ptr = ccgraph_get_component_row_ptr(&g, component);
    int *compo_col_ind = ccgraph_get_component_col_ind(&g, component);
    for (v_local=0; v_local<compo_nvertices; ++v_local) {
      for (i=compo_row_ptr[v_local]; i<compo_row_ptr[v_local+1]; ++i) {
        w_local = compo_col_ind[i];
        v_global = ccgraph_local_to_global(&g, component, v_local);
        w_global = ccgraph_local_to_global(&g, component, w_local);
        printf("  %d--%d\n", v_global, w_global);
      }
    }
    int compo_nedges = ccgraph_get_component_nedges(&g, component);
    printf("%d vertices\n", compo_nvertices);
    printf("%d edges\n", compo_nedges);
    printf("\n");
  }
  printf("\n");


  // Delete things from the connected components search.
  free(parent_ws);
  free(deck_ws);
  free(next_ws);
  free(component_labels);

  // Init the connected components search object.
  ccgraph_destroy(&g);

  // Delete the arrays that define the graph.
  free(row_ptr);
  free(col_ind);

  return 1;
}


int main()
{
  int nfails = 0;

  nfails += t0();

  if (nfails) {
    printf("failed testing: %d tests failed\n", nfails);
  } else {
    printf("success: all tests passed\n");
  }
}

