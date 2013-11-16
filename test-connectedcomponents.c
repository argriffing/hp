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

  // Init the connected components search object.
  CCGRAPH g;
  ccgraph_init(&g, nvertices, nedges);

  // Do the connected components search.
  ccgraph_compute(&g, row_ptr, col_ind, nvertices,
      parent_ws, deck_ws, next_ws,
      component_labels);

  // Delete things from the connected components search.
  free(parent_ws);
  free(deck_ws);
  free(next_ws);

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

