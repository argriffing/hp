#include "stdio.h"
#include "stdlib.h"

#include "graphgirth.h"


// Sparse matrix format conversion.
// List of lists to compressed row storage.
void lil_to_csr(int nvertices, int nedges,
    int *lil_in, int *row_ptr_out, int *col_ind_out)
{
  row_ptr_out[0] = 0;

  int v, e;
  int va, vb;
  int nfilled = 0;
  for (v=0; v<nvertices; ++v) {
    for (e=0; e<nedges; ++e) {
      int va = lil_in[e*2 + 0];
      int vb = lil_in[e*2 + 1];
      if (v == va) {
        col_ind_out[nfilled++] = vb;
      }
      if (v == vb) {
        col_ind_out[nfilled++] = va;
      }
    }
    row_ptr_out[v+1] = nfilled;
  }
}


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
  int *row_ptr = (int *) malloc((nvertices + 1) * sizeof(int));
  int *col_ind = (int *) malloc(2 * nedges * sizeof(int));
  lil_to_csr(nvertices, nedges, lil, row_ptr, col_ind);

  int root = 3;
  int *parent_ws = (int *) malloc(nvertices * sizeof(int));
  int *depth_ws = (int *) malloc(nvertices * sizeof(int));
  int *deck_ws = (int *) malloc(nvertices * sizeof(int));
  int *next_ws = (int *) malloc(nvertices * sizeof(int));
  int girth_ub = get_girth_ub(row_ptr, col_ind, nvertices, root,
      parent_ws, depth_ws, deck_ws, next_ws);
  free(parent_ws);
  free(depth_ws);
  free(deck_ws);
  free(next_ws);

  printf("root: %d\n", root);
  printf("girth upper bound: %d\n", girth_ub);

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

