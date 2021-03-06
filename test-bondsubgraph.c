#include "stdio.h"
#include "stdlib.h"
#include "assert.h"

#include "sparsetools.h"
#include "connectedcomponents.h"
#include "graphgirth.h"
#include "subsetsum.h"
#include "bondsubgraph.h"


// TODO copypasted from test-connectedcomponents
void _print_connected_component_graph(CCGRAPH *p)
{
  printf("connected component graphs:\n");
  int component, i;
  int v_local, v_global;
  int w_local, w_global;
  for (component=0; component<p->ncomponents; ++component) {
    printf("component %d:\n", component);
    int compo_nvertices = ccgraph_get_component_nvertices(p, component);
    int *compo_row_ptr = ccgraph_get_component_row_ptr(p, component);
    int *compo_col_ind = ccgraph_get_component_col_ind(p, component);
    for (v_local=0; v_local<compo_nvertices; ++v_local) {
      for (i=compo_row_ptr[v_local]; i<compo_row_ptr[v_local+1]; ++i) {
        w_local = compo_col_ind[i];
        v_global = ccgraph_local_to_global(p, component, v_local);
        w_global = ccgraph_local_to_global(p, component, w_local);
        printf("  %d--%d\n", v_global, w_global);
      }
    }
    int compo_nedges = ccgraph_get_component_nedges(p, component);
    printf("%d vertices\n", compo_nvertices);
    printf("%d edges\n", compo_nedges);
    printf("\n");
  }
  printf("\n");
}


int _solver_helper(int *lil, int nvertices, int nedges,
    int k, int *expected_solution)
{
  int failflag = 0;

  // Construct the csr representation of the graph.
  int *row_ptr = (int *) malloc((nvertices + 1) * sizeof(int));
  int *col_ind = (int *) malloc(2 * nedges * sizeof(int));
  lil_to_csr(nvertices, nedges, lil, row_ptr, col_ind);

  // Define some max values for initialization.
  int max_ncomponents = nvertices;
  int max_nvertices = nvertices;
  int max_target = nvertices;

  // Breadth first search, girth, and cycle allocation and initialization.
  BFS_WS bfs_ws;
  bfs_ws_init(&bfs_ws, nvertices);
  int *depth_ws = (int *) malloc(nvertices * sizeof(int));
  int *va_trace = (int *) malloc(nvertices * sizeof(int));
  int *vb_trace = (int *) malloc(nvertices * sizeof(int));

  int i, v;

  // Clear parent and depth arrays for girth and cycle related functions.
  for (v=0; v<nvertices; ++v) {
    bfs_ws.parent[v] = -1;
    depth_ws[v] = -1;
  }

  // Prepare to decompose the graph into connected components.
  int *component_labels = (int *) malloc(nvertices * sizeof(int));
  for (v=0; v<nvertices; ++v) {
    component_labels[v] = -1;
  }

  // Compute the connected components.
  CCGRAPH g;
  ccgraph_init(&g, nvertices, 2 * nedges);
  ccgraph_compute(&g, row_ptr, col_ind, nvertices, &bfs_ws, component_labels);

  // Initialize the subset sum solver workspaces.
  int *high = (int *) malloc(max_target * sizeof(int));
  int *low = (int *) malloc(max_target * sizeof(int));
  int *s3solution = (int *) malloc(max_ncomponents * sizeof(int));
  S3TABLE s3table;
  s3table_init(&s3table, max_ncomponents, max_target);
  s3table_clear(&s3table, max_ncomponents, max_target);

  // Init the densest k-subgraph solver.
  SOLVER solver;
  solver_init(&solver, max_nvertices, max_ncomponents);

  // Print the graph before trying to solve.
  printf("connected components before trying to solve densest k-subgraph:\n");
  _print_connected_component_graph(&g);

  // Get the densest k-subgraph.
  solver_toplevel(&solver,
      row_ptr, col_ind, &g, k,
      va_trace, vb_trace,
      &bfs_ws, depth_ws,
      &s3table, low, high, s3solution);

  // Check that the solver at least found the right number of vertices.
  assert(k == solver.nsolution);

  // Allocate some memory for comparing solutions.
  int *expected_set = (int *) calloc(max_nvertices, sizeof(int));
  int *observed_set = (int *) calloc(max_nvertices, sizeof(int));
  for (i=0; i<k; ++i) {
    v = expected_solution[i];
    expected_set[v] = 1;
  }
  for (i=0; i<solver.nsolution; ++i) {
    v = solver.solution[i];
    observed_set[v] = 1;
  }
  int mismatch_count = 0;
  for (v=0; v<max_nvertices; ++v) {
    printf("observed: %d  expected: %d\n",
        observed_set[v], expected_set[v]);
    if (observed_set[v] != expected_set[v]) {
      mismatch_count++;
    }
  }
  if (mismatch_count) {
    printf("%d mismatches\n", mismatch_count);
    failflag = 1;
  }

  // Explicitly count the number of directed edges in the induced subgraph.
  int j, w;
  int nedges_subgraph_directed = 0;
  for (i=0; i<solver.nsolution; ++i) {
    v = solver.solution[i];
    for (j=row_ptr[v]; j<row_ptr[v+1]; ++j) {
      w = col_ind[j];
      if (observed_set[v] && observed_set[w]) {
        nedges_subgraph_directed += 1;
      }
    }
  }
  assert(nedges_subgraph_directed % 2 == 0);
  int nedges_subgraph_undirected = nedges_subgraph_directed / 2;
  if (nedges_subgraph_undirected != solver.score) {
    printf("failed to reproduce the reported score:\n");
    printf("reported score: %d\n", solver.score);
    printf("computed score: %d\n", nedges_subgraph_undirected);
    failflag = 1;
  }

  // Print the graph after some re-ordering.
  // Note that this does not reflect the re-ordering
  // of components within the solver.
  printf("connected components after having tried to solve:\n");
  _print_connected_component_graph(&g);

  // Report the best solution.
  printf("computed densest k-subgraph vertex set:\n");
  for (v=0; v<solver.nsolution; ++v) {
    printf("%d ", solver.solution[v]);
  }
  printf("\n");

  return failflag;
}



int _girth_conn_helper(int *lil, int nvertices, int nedges,
    int girth, int *expected_small_cycle)
{
  int i, v;
  int failflag = 0;

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

  // Initialize the mask for the expected small cycle.
  int *small_cycle_mask = (int *) calloc(nvertices, sizeof(int));
  for (i=0; i<girth; ++i) {
    v = expected_small_cycle[i];
    small_cycle_mask[v] = 1;
  }

  // Initialize component labels.
  int *component_labels = (int *) malloc(nvertices * sizeof(int));
  for (v=0; v<nvertices; ++v) {
    component_labels[v] = -1;
  }

  // Compute the connected components.
  CCGRAPH g;
  ccgraph_init(&g, nvertices, 2 * nedges);
  ccgraph_compute(&g, row_ptr, col_ind, nvertices,
      &bfs_ws, component_labels);

  if (g.ncomponents != 1) {
    printf("incorrect number of connected components\n");
    printf("expected: %d\n", 1);
    printf("computed: %d\n", g.ncomponents);
    printf("\n");
    failflag = 1;
  }

  printf("graph before moving the small cycle to the front:\n");
  _print_connected_component_graph(&g);

  // Try to move the smallest cycle to the front.
  int component = 0;
  int computed_girth =  _move_smallest_cycle_to_front(
      row_ptr, col_ind,
      &g, component,
      va_trace, vb_trace,
      &bfs_ws, depth_ws);

  if (girth != computed_girth) {
    printf("incorrect computed girth:\n");
    printf("expected: %d\n", girth);
    printf("computed: %d\n", computed_girth);
    printf("\n");
    failflag = 1;
  }

  printf("graph after moving the small cycle to the front:\n");
  _print_connected_component_graph(&g);

  // Free some memory.
  ccgraph_destroy(&g);
  bfs_ws_destroy(&bfs_ws);
  free(component_labels);
  free(va_trace);
  free(vb_trace);
  free(small_cycle_mask);
  free(depth_ws);
  free(row_ptr);
  free(col_ind);

  return failflag;
}




// TODO this is copypasted from the girth test.
//
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

int test_three_cycles()
{
  // Check if we can solve the subset sum problem for this graph.
  // There is only one solution.
  // 
  // 0 -- 1   7 -- 6 -- 5 -- 4
  // |    |   |    |
  // 3 -- 2   8   13
  //          |    |
  //          9   12
  //          |    |
  //         10 - 11
  //
  // 14 -- 15 -- 16
  //  |           |
  // 17 -- 18 -- 19
  //
  int lil[] = {
    0, 1,
    1, 2,
    2, 3,
    3, 0,

    4, 5,
    5, 6,
    6, 7,
    7, 8,
    8, 9,
    9, 10,
    10, 11,
    11, 12,
    12, 13,
    13, 6,
  
    14, 15,
    15, 16,
    16, 19,
    19, 18,
    18, 17,
    17, 14,
  };
  int nvertices = 20;
  int nedges = 20;
  int k = 15;
  int expected_solution[] = {
    14, 15, 16, 17, 18, 19,
    5, 6, 7, 8, 9, 10, 11, 12, 13};
  return _solver_helper(lil, nvertices, nedges, k, expected_solution);
}


int test_menagerie()
{
  //
  //     20 -- 22   19 -- 23 -- 21
  //                       |
  //                      17
  //     0 --- 1
  //      \   /
  //       \ /       13 -- 14 -- 15     9
  //        2               |     |     |
  //                       16 -- 12 -- 10 -- 11
  // 
  //  6 -- 7 -- 8
  //  |    |    |       18
  //  3 -- 4 -- 5
  //
  //
  int lil[] = {
    3, 4,
    4, 5,
    6, 7,
    7, 8,
    6, 3,
    7, 4,
    8, 5,

    0, 1,
    1, 2,
    2, 0,

    20, 22,

    19, 23,
    23, 21,
    23, 17,

    13, 14,
    14, 15,
    14, 16,
    16, 12,
    15, 12,
    11, 10,
    10, 9,
    10, 12,
  };
  int nvertices = 24;
  int nedges = 22;
  int k = 21;
  int expected_solution[] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
    13, 14, 15, 16, 17, 19, 21, 23};
  return _solver_helper(lil, nvertices, nedges, k, expected_solution);
}

int test_duelling_cycles()
{
  // 0 -- 1   3   6 -- 9
  // |   /   /|   |    |
  // |  /   / |   7 -- 10
  // | /   /  |   |    |
  // |/   /   |   8 -- 11
  // 2   4 -- 5
  //
  int lil[] = {
    0, 1,
    1, 2,
    2, 0,

    3, 4,
    4, 5,
    5, 3,

    6, 9,
    7, 10,
    8, 11,
    6, 7,
    7, 8,
    9, 10,
    10, 11,
  };
  int nvertices = 12;
  int nedges = 13;
  int k = 6;
  int expected_solution[] = {6, 7, 8, 9, 10, 11};
  return _solver_helper(lil, nvertices, nedges, k, expected_solution);
}

int test_pure_isolation()
{
  // 0           4
  //      1
  //  2      3
  //
  int *lil;
  int nvertices = 5;
  int nedges = 0;
  int k = 5;
  int expected_solution[] = {0, 1, 2, 3, 4};
  return _solver_helper(lil, nvertices, nedges, k, expected_solution);
}

// This test has a huge component that gives us indigestion.
int test_indigestion()
{
  //     4       17
  //  5                   19           9
  //        20                  10
  //    22        12                          14
  //                        15
  //   3 -- 6
  //                               2 ----- 23
  //                                \
  //                                 \
  //                                  \
  //                                   7
  //
  //  16 -- 25 -- 18 -- 21 -- 24 -- 1
  //   |    |                       |
  //  13 -- 0 --- 11 -------------- 8
  //
  int lil[] = {
    23, 2,
    2, 7,
    3, 6,
    16, 25,
    25, 18,
    18, 21,
    21, 24,
    24, 1,
    1, 8,
    8, 11,
    11, 0,
    0, 13,
    13, 16,
    0, 25,
  };
  int nvertices = 26;
  int nedges = 14;
  int k = 6;
  int expected_solution[] = {0, 13, 16, 25, 11, 18};
  return _solver_helper(lil, nvertices, nedges, k, expected_solution);
}


int main()
{
  int nfails = 0;

  nfails += test_small_cycle_1();
  nfails += test_three_cycles();
  nfails += test_menagerie();
  nfails += test_duelling_cycles();
  nfails += test_pure_isolation();
  nfails += test_indigestion();

  if (nfails) {
    printf("failed testing: %d tests failed\n", nfails);
  } else {
    printf("success: all tests passed\n");
  }
}

