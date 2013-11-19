#include "stdio.h"
#include "stdlib.h"

#include "sparsetools.h"
#include "graphgirth.h"

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

  // Allocate csr memory and convert to csr format.
  int *row_ptr = (int *) malloc((nvertices + 1) * sizeof(int));
  int *col_ind = (int *) malloc(2 * nedges * sizeof(int));
  lil_to_csr(nvertices, nedges, lil, row_ptr, col_ind);

  // Allocate memory for the search for the smallest cycle.
  BFS_WS bfs_ws;
  bfs_ws_init(&bfs_ws, nvertices);
  int *depth_ws = (int *) malloc(nvertices * sizeof(int));

  // Compute girth with and without the connectivity assumption.
  int computed_girth_conn = get_girth_conn(row_ptr, col_ind, nvertices,
      &bfs_ws, depth_ws);
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

  // Free some memory.
  bfs_ws_destroy(&bfs_ws);
  free(depth_ws);
  free(row_ptr);
  free(col_ind);

  return failflag;
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

  // Allocate csr memory and convert to csr format.
  int *row_ptr = (int *) malloc((nvertices + 1) * sizeof(int));
  int *col_ind = (int *) malloc(2 * nedges * sizeof(int));
  lil_to_csr(nvertices, nedges, lil, row_ptr, col_ind);

  // Allocate memory for the search for the smallest cycle.
  BFS_WS bfs_ws;
  bfs_ws_init(&bfs_ws, nvertices);
  int *depth_ws = (int *) malloc(nvertices * sizeof(int));

  // Compute girth with and without the connectivity assumption.
  int computed_girth_conn = get_girth_conn(row_ptr, col_ind, nvertices,
      &bfs_ws, depth_ws);
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

  // Free some memory.
  bfs_ws_destroy(&bfs_ws);
  free(depth_ws);
  free(row_ptr);
  free(col_ind);

  return failflag;
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

  // Allocate csr memory and convert to csr format.
  int *row_ptr = (int *) malloc((nvertices + 1) * sizeof(int));
  int *col_ind = (int *) malloc(2 * nedges * sizeof(int));
  lil_to_csr(nvertices, nedges, lil, row_ptr, col_ind);

  // Allocate memory for the search for the smallest cycle.
  BFS_WS bfs_ws;
  bfs_ws_init(&bfs_ws, nvertices);
  int *depth_ws = (int *) malloc(nvertices * sizeof(int));

  // Smoke test for girth upper bound.
  int root = 3;
  int computed_girth_ub = get_girth_ub(row_ptr, col_ind, nvertices, root,
      &bfs_ws, depth_ws);

  // Compute girth with and without the connectivity assumption.
  int computed_girth_conn = get_girth_conn(row_ptr, col_ind, nvertices,
      &bfs_ws, depth_ws);
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

  // Free some memory.
  bfs_ws_destroy(&bfs_ws);
  free(depth_ws);
  free(row_ptr);
  free(col_ind);

  return failflag;
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

