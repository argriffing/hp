#include "stdio.h"
#include "stdlib.h"

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

  // Declare vertex for iteration.
  int v;

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

  // Report the best solution.
  printf("computed densest k-subgraph vertex set:\n");
  for (v=0; v<solver.nsolution; ++v) {
    printf("%d ", solver.solution[v]);
  }
  printf("\n");

  return 1;
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

int main()
{
  int nfails = 0;

  nfails += test_small_cycle_1();
  nfails += test_three_cycles();

  if (nfails) {
    printf("failed testing: %d tests failed\n", nfails);
  } else {
    printf("success: all tests passed\n");
  }
}

