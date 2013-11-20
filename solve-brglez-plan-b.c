#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "string.h"

#include "sparsetools.h"
#include "connectedcomponents.h"
#include "graphgirth.h"
#include "subsetsum.h"
#include "bondsubgraph.h"


int solve_potential_bond_graph(
    const int *row_ptr, const int *col_ind, int nvertices, int k)
{
  int failflag = 0;

  int max_nvertices = nvertices;
  int max_ncomponents = nvertices;
  int max_target = nvertices;
  int max_nedges_undirected = nvertices + 1;
  int max_nedges_directed = max_nedges_undirected * 2;

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
  ccgraph_init(&g, max_nvertices, max_nedges_directed);
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

  // Get the densest k-subgraph.
  solver_toplevel(&solver,
      row_ptr, col_ind, &g, k,
      va_trace, vb_trace,
      &bfs_ws, depth_ws,
      &s3table, low, high, s3solution);

  int *binary_solution = (int *) calloc(nvertices, sizeof(int));
  for (i=0; i<solver.nsolution; ++i) {
    v = solver.solution[i];
    binary_solution[v] = 1;
  }

  // Print the score and the binary representation.
  printf("%d\n", solver.score);
  for (v=0; v<nvertices; ++v) {
    printf("%d", binary_solution[v]);
  }
  printf("\n");

  // Free the memory.
  bfs_ws_destroy(&bfs_ws);
  ccgraph_destroy(&g);
  s3table_destroy(&s3table);
  solver_destroy(&solver);
  free(binary_solution);
  free(high);
  free(low);
  free(s3solution);
  free(component_labels);
  free(va_trace);
  free(vb_trace);
  free(depth_ws);

  return failflag;
}


// dlluruluurdrdrurddl
// rdruruluuurululdldrdddl
// lurulurulldlulddrdddruuu
int main(int argc, char *argv[])
{
  if (argc != 3) {
    printf("Usage:\n");
    printf("%s 10 dlluruluurdrdrurddl\n", argv[0]);
    return 1;
  }
  int i;
  char c;

  int k = atoi(argv[1]);
  int path_length = strlen(argv[2]);

  // Check the input before we get too far.

  if (k < 2) {
    printf("Error: the sequence weight should be at least 2\n");
    return 1;
  }

  if (path_length < 3) {
    printf("Error: the walk should take at least 3 steps\n");
    return 1;
  }

  for (i=0; i<path_length; ++i) {
    c = argv[2][i];
    if (c != 'l' && c != 'r' && c != 'u' && c != 'd') {
      printf("Error: each step of the walk must be one of {l, r, u, d}\n");
      return 1;
    }
  }

  // Define the number of vertices in the potential bond graph,
  // and an upper bound on the number of edges in the potential bond graph.
  int nvertices = path_length + 1;
  int nedges_max_undirected = nvertices + 1;
  int nedges_max_directed = 2 * nedges_max_undirected;

  // Allocate room for the csr graph.
  int *row_ptr = (int *) malloc((nvertices + 1) * sizeof(int));
  int *col_ind = (int *) malloc(nedges_max_directed * sizeof(int));

  // Allocate room for x and y coordinates per vertex.
  int *x = (int *) malloc(nvertices * sizeof(int));
  int *y = (int *) malloc(nvertices * sizeof(int));

  // Initialize the per-vertex coordinates.
  x[0] = 0;
  y[0] = 0;
  int step;
  for (step=0; step<path_length; ++step) {
    x[step + 1] = x[step];
    y[step + 1] = y[step];
    c = argv[2][step];
    if (c == 'l') x[step + 1]--;
    if (c == 'r') x[step + 1]++;
    if (c == 'u') y[step + 1]--;
    if (c == 'd') y[step + 1]++;
  }

  // Initialize the csr graph of potential bonds.
  row_ptr[0] = 0;
  int j;
  for (i=0; i<nvertices; ++i) {
    row_ptr[i+1] = row_ptr[i];
    for (j=0; j<nvertices; ++j) {
      if (abs(j-i) > 1) {
        if (abs(x[j] - x[i]) + abs(y[j] - y[i]) == 1) {
          col_ind[row_ptr[i+1]++] = j;
        }
      }
    }
  }

  // Free the coordinates per vertex.
  free(x);
  free(y);

  // Solve a densest k-subgraph problem on the potential bond graph.
  int failflag = solve_potential_bond_graph(row_ptr, col_ind, nvertices, k);

  // Free the csr graph.
  free(row_ptr);
  free(col_ind);

  // More or less success.
  return failflag;
}

