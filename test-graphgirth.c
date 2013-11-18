#include "stdio.h"
#include "stdlib.h"

#include "sparsetools.h"
#include "graphgirth.h"

int t0()
{
  // This graph is connected.
  // It is like a theta graph with a tree attached onto it.
  // The smallest cycle of this graph is {7, 8, 10, 9}.
  // If you root at early nodes like < 5 then you will
  // not find that cycle using a topo sort detection.
  int lil[] = {
    0, 1,
    1, 2,
    1, 3,
    1, 4,
    4, 5,
    4, 11,
    5, 6,
    6, 7,
    7, 9,
    7, 8,
    8, 10,
    9, 10,
    11, 12,
    12, 10,
  };
  int nvertices = 13;
  int nedges = 14;

  // Allocate csr memory and convert to csr format.
  int *row_ptr = (int *) malloc((nvertices + 1) * sizeof(int));
  int *col_ind = (int *) malloc(2 * nedges * sizeof(int));
  lil_to_csr(nvertices, nedges, lil, row_ptr, col_ind);

  int root = 3;

  // Allocate memory for the search for the smallest cycle.
  BFS_WS bfs_ws;
  bfs_ws_init(&bfs_ws, nvertices);
  int *depth_ws = (int *) malloc(nvertices * sizeof(int));

  // Compute girth.
  int girth_ub = get_girth_ub(row_ptr, col_ind, nvertices, root,
      &bfs_ws, depth_ws);
  int girth = get_girth(row_ptr, col_ind, nvertices,
      &bfs_ws, depth_ws);

  printf("root: %d\n", root);
  printf("girth upper bound: %d\n", girth_ub);
  printf("\n");

  printf("girth: %d\n", girth);
  printf("\n");
  
  // Free girth calculation memory.
  bfs_ws_destroy(&bfs_ws);
  free(depth_ws);

  // Free csr memory.
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

