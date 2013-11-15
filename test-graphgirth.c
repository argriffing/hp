#include "graphgirth.h"

// Sparse matrix format conversion.
// List of lists to compressed row storage.
void lil_to_csr(int nvertices, int nedges,
    int *lil_in, int *row_ptr_out, int *col_ind_out)
{
}


int t0()
{
  // This graph is connected.
  // It is like a theta graph with a tree attached onto it.
  // The smallest cycle of this graph is {7, 8, 10, 9}.
  // If you root at early nodes like < 5 then you will
  // not find that cycle using a topo sort detection.
  int lil = {
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
  }
  int nvertices = 13;
  int nedges = 14;
  int *row_ptr = (int *) malloc((nvertices + 1) * sizeof(int));
  int *col_ind = (int *) malloc(nedges * sizeof(int));
  lil_to_csr(nvertices, nedges, lil, row_ptr, col_ind);

  free(row_ptr);
  free(col_ind);
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

