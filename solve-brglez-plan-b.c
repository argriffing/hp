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
    const int *row_ptr, const int *col_ind, int nvertices, int k,
    int *binary_solution, int *pbest_score)
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

  // Get the number of connected components that have cycles.
  // The purpose is to limit the size of the dynamic programming table
  // that is preallocated for solving the subset sum problem.
  int cyclic_component_count = 0;
  int c;
  for (c=0; c<g.ncomponents; ++c)
  {
    int component_nedges_undirected = ccgraph_get_component_nedges_undirected(
        &g, c);
    int component_nvertices = ccgraph_get_component_nvertices(
        &g, c);
    if (component_nedges_undirected >= component_nvertices) {
      cyclic_component_count++;
    }
  }

  // Initialize the subset sum solver workspaces.
  int *high = (int *) malloc(max_target * sizeof(int));
  int *low = (int *) malloc(max_target * sizeof(int));
  int *s3solution = (int *) malloc(cyclic_component_count * sizeof(int));
  S3TABLE s3table;
  s3table_init(&s3table, cyclic_component_count, max_target);
  s3table_clear(&s3table, cyclic_component_count, max_target);

  // Init the densest k-subgraph solver.
  SOLVER solver;
  solver_init(&solver, max_nvertices, max_ncomponents);

  // Get the densest k-subgraph.
  solver_toplevel(&solver,
      row_ptr, col_ind, &g, k,
      va_trace, vb_trace,
      &bfs_ws, depth_ws,
      &s3table, low, high, s3solution);

  // Fill the output values.
  for (i=0; i<solver.nsolution; ++i) {
    v = solver.solution[i];
    binary_solution[v] = 1;
  }
  *pbest_score = solver.score;

  // Free the memory.
  bfs_ws_destroy(&bfs_ws);
  ccgraph_destroy(&g);
  s3table_destroy(&s3table);
  solver_destroy(&solver);
  free(high);
  free(low);
  free(s3solution);
  free(component_labels);
  free(va_trace);
  free(vb_trace);
  free(depth_ws);

  return failflag;
}

// TODO: use hashing to construct the potential bond graph in O(N)

// Construct the potential bond graph using an explicit bounding box.
// Depending on the efficiency of the box,
// the complexity is between O(N) and O(N^2) time and space.
int construct_bond_graph_using_box(
    const int *x, const int *y, int nvertices,
    int *row_ptr, int *col_ind)
{
  assert(nvertices);
  int i;

  // Find the bounds of the box.
  int xmin = x[0];
  int xmax = x[0];
  int ymin = y[0];
  int ymax = y[0];
  for (i=0; i<nvertices; ++i) {
    if (x[i] > xmax) xmax = x[i];
    if (x[i] < xmin) xmin = x[i];
    if (y[i] > ymax) ymax = y[i];
    if (y[i] < ymin) ymin = y[i];
  }

  // Allocate a bordered table.
  int width = (xmax - xmin + 1) + 2;
  int height = (ymax - ymin + 1) + 2;
  int area = width * height;
  int xbase = xmin - 1;
  int ybase = ymin - 1;
  int *table = (int *) malloc(area * sizeof(int));
  if (!table) {
    printf("Error: memory error allocating bounding box\n");
    return 1;
  }
  memset(table, -1, area * sizeof(int));

  // Fill the table, checking for collisions.
  int idx;
  int ncollisions = 0;
  for (i=0; i<nvertices; ++i) {
    idx = (y[i] - ybase)*width + (x[i] - xbase);
    if (table[idx] != -1) {
      printf("Error: the conformation includes self-intersections\n");
      free(table);
      return 1;
    }
    table[idx] = i;
  }

  // Add the potential bond edges.
  int dx[] = {1, 0, -1, 0};
  int dy[] = {0, -1, 0, 1};
  int j;
  int jdx;
  row_ptr[0] = 0;
  for (i=0; i<nvertices; ++i) {
    row_ptr[i+1] = row_ptr[i];
    idx = (y[i] - ybase)*width + (x[i] - xbase);
    for (j=0; j<4; ++j) {
      jdx = (y[i] + dy[j] - ybase)*width + (x[i] + dx[j] - xbase);
      if (table[jdx] >= 0 && abs(table[jdx] - table[idx]) > 1) {
        col_ind[row_ptr[i+1]++] = table[jdx];
      }
    }
  }

  // Free the table.
  free(table);
  return 0;
}


// Construct the potential bond graph in O(N^2) time.
int construct_bond_graph(
    const int *x, const int *y, int nvertices,
    int *row_ptr, int *col_ind)
{
  int i, j;
  int ncollisions = 0;
  row_ptr[0] = 0;
  for (i=0; i<nvertices; ++i) {
    row_ptr[i+1] = row_ptr[i];
    for (j=0; j<nvertices; ++j) {
      if (i != j && x[i] == x[j] && y[i] == y[j]) {
        ncollisions++;
      }
      if (abs(j-i) > 1) {
        if (abs(x[j] - x[i]) + abs(y[j] - y[i]) == 1) {
          col_ind[row_ptr[i+1]++] = j;
        }
      }
    }
  }
  if (ncollisions) {
    printf("Error: the conformation includes self-intersections\n");
  }
  return ncollisions;
}


// dlluruluurdrdrurddl
// rdruruluuurululdldrdddl
// lurulurulldlulddrdddruuu
int main(int argc, char *argv[])
{
  if (argc != 2) {
    printf("Usage:\n");
    printf("%s 10\n", argv[0]);
    return 1;
  }
  int i;
  char c;

  int k = atoi(argv[1]);
  if (k < 2) {
    printf("Error: the sequence weight should be at least 2\n");
    return 1;
  }

  // Allocate room for the walk to be read from stdin.
  int max_walk_len = 1000000;
  int max_walk_size = max_walk_len + 1;
  char *walk = (char *) calloc(max_walk_size, sizeof(char));

  // Read the walk from stdin.
  fgets(walk, max_walk_size, stdin);
  int raw_walk_length = strlen(walk);

  // Check that the input is not too big.
  if (raw_walk_length == max_walk_len) {
    printf("Error: please use a walk input less than %d bytes\n", max_walk_len);
    return 1;
  }

  // Remove whitespace.
  char *pfront = walk;
  char *pback = walk;
  while (*pback) {
    c = *(pback++);
    if (!isspace(c)) {
      *(pfront++) = c;
    }
  }
  *pfront = 0;
  int walk_length = strlen(walk);

  // Check the input before we get too far.

  if (walk_length < 3) {
    printf("Error: the walk should take at least 3 steps\n");
    return 1;
  }

  for (i=0; i<walk_length; ++i) {
    c = walk[i];
    if (c != 'l' && c != 'r' && c != 'u' && c != 'd') {
      printf("Error: each step of the walk must be one of {l, r, u, d}\n");
      return 1;
    }
  }

  // Define the number of vertices in the potential bond graph,
  // and an upper bound on the number of edges in the potential bond graph.
  int nvertices = walk_length + 1;
  int nedges_max_undirected = nvertices + 1;
  int nedges_max_directed = 2 * nedges_max_undirected;

  // If the requested weight is larger than the number of vertices,
  // then this is a problem.
  if (k > nvertices) {
    printf("Error: the requested sequence weight ");
    printf("is greater than the sequence length\n");
    free(walk);
    return 1;
  }

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
  for (step=0; step<walk_length; ++step) {
    x[step + 1] = x[step];
    y[step + 1] = y[step];
    c = walk[step];
    if (c == 'l') x[step + 1]--;
    if (c == 'r') x[step + 1]++;
    if (c == 'u') y[step + 1]--;
    if (c == 'd') y[step + 1]++;
  }

  // Free the walk memory.
  free(walk);

  // Initialize the csr graph of potential bonds.
  int construction_failure = construct_bond_graph_using_box(
      x, y, nvertices, row_ptr, col_ind);

  // Free the coordinates per vertex.
  free(x);
  free(y);

  if (construction_failure) {
    free(row_ptr);
    free(col_ind);
    return construction_failure;
  }

  // Solve a densest k-subgraph problem on the potential bond graph.
  int binary_solution[nvertices];
  memset(binary_solution, 0, sizeof binary_solution);

  int best_score = -1;
  int failflag = solve_potential_bond_graph(row_ptr, col_ind, nvertices, k,
      binary_solution, &best_score);

  // Print the score and the binary representation.
  int v;
  printf("%d\n", best_score);
  for (v=0; v<nvertices; ++v) {
    printf("%d", binary_solution[v]);
  }
  printf("\n");

  // Free the csr graph.
  free(row_ptr);
  free(col_ind);

  // More or less success.
  return failflag;
}

