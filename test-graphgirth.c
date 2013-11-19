#include "stdio.h"
#include "stdlib.h"

#include "sparsetools.h"
#include "graphgirth.h"


int _girth_conn_helper(int *lil, int nvertices, int nedges,
    int girth, int *expected_small_cycle)
{
  int i, v;

  // Allocate csr memory and convert to csr format.
  int *row_ptr = (int *) malloc((nvertices + 1) * sizeof(int));
  int *col_ind = (int *) malloc(2 * nedges * sizeof(int));
  lil_to_csr(nvertices, nedges, lil, row_ptr, col_ind);

  // Allocate memory for the search for the smallest cycle.
  BFS_WS bfs_ws;
  bfs_ws_init(&bfs_ws, nvertices);
  int *depth_ws = (int *) malloc(nvertices * sizeof(int));

  // Initialize parent and depth to -1 for all vertices.
  for (v=0; v<nvertices; ++v) {
    bfs_ws.parent[v] = -1;
    depth_ws[v] = -1;
  }

  // Allocate memory to hold the smallest cycle.
  int *va_trace = (int *) malloc(nvertices * sizeof(int));
  int *vb_trace = (int *) malloc(nvertices * sizeof(int));
  int *small_cycle = (int *) malloc(nvertices * sizeof(int));

  // Initialize the mask for the expected small cycle.
  int *small_cycle_mask = (int *) calloc(nvertices, sizeof(int));
  for (i=0; i<girth; ++i) {
    v = expected_small_cycle[i];
    small_cycle_mask[v] = 1;
  }

  // Compute girth with and without the connectivity assumption.
  int computed_girth_conn = -1;
  int girth_witness = -1;
  get_girth_and_vertex_conn(row_ptr, col_ind, nvertices,
      &bfs_ws, depth_ws,
      &computed_girth_conn, &girth_witness);
  int computed_girth = get_girth(row_ptr, col_ind, nvertices,
      &bfs_ws, depth_ws);

  int failflag = 0;

  // Check the girth, assuming that the graph is connected.
  if (computed_girth_conn != girth) {
    printf("failed to compute the girth correctly,\n");
    printf("under the assumption of graph connnectivity\n");
    printf("expected girth: %d\n", girth);
    printf("computed girth: %d\n", computed_girth_conn);
    printf("\n");
    failflag = 1;
  }

  // Check the girth.
  if (computed_girth != girth) {
    printf("failed to compute the girth correctly\n");
    printf("expected girth: %d\n", girth);
    printf("computed girth: %d\n", computed_girth);
    printf("\n");
    failflag = 1;
  }

  // Extract the smallest cycle using the witness vertex
  // of the girth calculation.
  int small_cycle_length = -1;
  if (girth_witness >= 0) {
    get_smallest_cycle_ub(
        row_ptr, col_ind, girth_witness,
        &bfs_ws, depth_ws,
        va_trace, vb_trace,
        small_cycle, &small_cycle_length);
  }

  // Check the extraction of the small cycle.
  if (small_cycle_length != girth) {
    printf("failed to extract a cycle with the correct girth\n");
    printf("expected girth: %d\n", girth);
    printf("computed girth: %d\n", small_cycle_length);
    printf("\n");
    failflag = 1;
  }

  // If the length of the extracted small cycle is correct,
  // then check that the extracted cycle is correct.
  if (girth_witness >= 0 && small_cycle_length == girth) {
    int small_cycle_fail = 0;
    for (i=0; i<small_cycle_length; ++i) {
      v = small_cycle[i];
      if (small_cycle_mask[v] == 0) {
        small_cycle_fail = 1;
        failflag = 1;
      }
    }
    if (small_cycle_fail) {
      printf("incorrect small cycle vertex set\n");
      printf("expected: ");
      for (i=0; i<nvertices; ++i) {
        if (small_cycle_mask[i]) {
          printf("%d ", i);
        }
      }
      printf("\n");
      printf("computed: ");
      for (i=0; i<small_cycle_length; ++i) {
        v = small_cycle[i];
        printf("%d ", v);
      }
      printf("\n");
    }
  }

  // Free some memory.
  bfs_ws_destroy(&bfs_ws);
  free(va_trace);
  free(vb_trace);
  free(small_cycle);
  free(small_cycle_mask);
  free(depth_ws);
  free(row_ptr);
  free(col_ind);

  return failflag;
}


int test_small_cycle_1()
{
  // Check if we can find the small cycle
  // in a graph with multiple overlapping cycles.
  //
  //   0 ------- 1 ------- 2 ------- 3
  //   |                             |
  //   4 -- 5 -- 6 -- 7 -- 8 -- 9 -- 10
  //   |                             |
  //  11 ------ 12 ------- 13 ------ 14
  //
  int lil[] = {
    0, 1,
    1, 2,
    2, 3,
    4, 5,
    5, 6,
    6, 7,
    7, 8,
    8, 9,
    9, 10,
    11, 12,
    12, 13,
    13, 14,
    0, 4,
    4, 11,
    3, 10,
    10, 14,
  };
  int nvertices = 15;
  int nedges = 16;
  int girth = 10;
  int small_cycle[] = {0, 1, 2, 3, 11, 12, 13, 14, 4, 10};
  return _girth_conn_helper(lil, nvertices, nedges, girth, small_cycle);
}

int test_small_cycle_0()
{
  // Check if we can find the small cycle in a hollow-barbell graph.
  //
  //   0--1
  //   |  |
  //   3--2
  //   |
  //   4
  //   |
  //   5
  //   |
  //   6--7
  //   |  |
  //  11  8
  //   |  |
  //  10--9
  //
  int lil[] = {
    0, 1,
    1, 2,
    2, 3,
    3, 0,
    3, 4,
    4, 5,
    5, 6,
    6, 7,
    7, 8,
    8, 9,
    9, 10,
    10, 11,
    11, 6,
  };
  int nvertices = 12;
  int nedges = 13;
  int girth = 4;
  int small_cycle[] = {0, 1, 2, 3};
  return _girth_conn_helper(lil, nvertices, nedges, girth, small_cycle);
}




int t0()
{
  // This graph is connected.
  // It is like a theta graph with a tree attached onto it.
  // The smallest cycle of this graph is {7, 8, 10, 9},
  // so the girth is 4.
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
  int girth = 4;
  int small_cycle[] = {7, 8, 9, 10};
  return _girth_conn_helper(lil, nvertices, nedges, girth, small_cycle);
}


int main()
{
  int nfails = 0;

  nfails += t0();
  nfails += test_small_cycle_0();
  nfails += test_small_cycle_1();

  if (nfails) {
    printf("failed testing: %d tests failed\n", nfails);
  } else {
    printf("success: all tests passed\n");
  }
}

