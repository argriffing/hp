#include "assert.h"
#include "stdio.h"

#include "breadthfirst.h"
#include "graphgirth.h"


// Get the two nodes furthest away from vertex r
// in the smallest cycle detected in a topo sort rooted at r.
// The breadth first search workspace has .parent, .deck, and .next.
void _girth_ub_helper(
    const int *row_ptr, const int *col_ind, int nvertices, int r,
    BFS_WS *bfs_ws, int *depth_ws,
    int *min_length, int *va_out, int *vb_out
    )
{
  int i, v;

  // The min cycle length is unknown,
  // and the nodes at the far end of the cycle are unknown.
  *min_length = -1;
  *va_out = -1;
  *vb_out = -1;

  // No vertices have been visited yet.
  // A parent of -2 means root.
  // A parent of -1 means the vertex has not yet been visited.
  // A parent >= 0 is the parent vertex index.
  // A depth of -1 means the vertex has not yet been visited.
  // A depth of >= 0 is the distance from the root.
  for (v=0; v<nvertices; ++v) {
    bfs_ws->parent[v] = -1;
    depth_ws[v] = -1;
  }
  bfs_ws->parent[r] = -2;
  depth_ws[r] = 0;

  // The root vertex is on deck for breadth first search.
  // The root has depth zero.
  bfs_ws->deck[0] = r;
  int ndeck = 1;

  // Stop after the first breadth first search iteration that finds a cycle.
  // Note that we cannot stop immediately after the first cycle is found,
  // because cycles found in the same iteration may differ in length
  // by up to one unit.
  int depth = 0;
  while (ndeck && (*min_length < 0))
  {
    int nnext = 0;
    int ideck;
    for (ideck=0; ideck<ndeck; ++ideck) {
      v = bfs_ws->deck[ideck];
      for (i=row_ptr[v]; i<row_ptr[v+1]; ++i) {
        int w = col_ind[i];
        if (w != bfs_ws->parent[v]) {

          // Look at each vertex one step outside of the current shell.
          // If the neighbor has already been spotted,
          // then check if the cycle length is the smallest known.
          // If the neighbor has not been spotted
          // then set its parent and depth and add it to the shell.
          if (bfs_ws->parent[w] >= 0) {
            int cycle_length = depth_ws[v] + depth_ws[w] + 1;
            if (*min_length < 0 || cycle_length < *min_length) {
              //printf("(v, w): (%d, %d)\n", v, w);
              //printf("depth of v: %d\n", depth_ws[v]);
              //printf("depth of w: %d\n", depth_ws[w]);
              //printf("cycle length %d\n", cycle_length);
              //printf("\n");
              *min_length = cycle_length;
              *va_out = v;
              *vb_out = w;
            }
          } else {
            bfs_ws->parent[w] = v;
            depth_ws[w] = depth_ws[v] + 1;
            bfs_ws->next[nnext++] = w;
          }
        }
      }
    }

    // Put the next array on deck.
    int *tmp = bfs_ws->deck;
    bfs_ws->deck = bfs_ws->next;
    bfs_ws->next = tmp;
    ndeck = nnext;
  }
}

// Get the length of the smallest cycle detected from topo sort rooted at r.
// This is an upper bound on the girth of the graph.
// The _ws suffixed arrays should be of length nvertices
// and are for temporary storage.
int get_girth_ub(
    const int *row_ptr, const int *col_ind, int nvertices, int r,
    BFS_WS *bfs_ws, int *depth_ws
    )
{
  int min_length, va, vb;
  _girth_ub_helper(row_ptr, col_ind, nvertices, r,
      bfs_ws, depth_ws,
      &min_length, &va, &vb);
  return min_length;
}

// Get the smallest cycle detected from topo sort rooted at r.
// The length of this cycle is an upper bound
// on the smallest cycle in the graph.
int get_smallest_cycle_ub(
    const int *row_ptr, const int *col_ind, int nvertices, int r,
    BFS_WS *bfs_ws, int *depth_ws,
    int *cycle_out, int *ncycle_out
    )
{
  int min_length, va, vb;
  _girth_ub_helper(row_ptr, col_ind, nvertices, r,
      bfs_ws, depth_ws, 
      &min_length, &va, &vb);

  // Reuse the deck_ws and next_ws buffers for other purposes.
  int *va_trace = bfs_ws->deck;
  int *vb_trace = bfs_ws->next;

  // Use the parent_ws output from the helper function
  // to trace both vertices back to the root.
  int i, v;
  int na, nb;

  v = va;
  na = 0;
  while (v >= 0) {
    va_trace[na++] = v;
    v = bfs_ws->parent[v];
  }

  v = vb;
  nb = 0;
  while (v >= 0) {
    vb_trace[nb++] = v;
    v = bfs_ws->parent[v];
  }

  // Check that both traces include at least two vertices.
  assert(na >= 2);
  assert(nb >= 2);

  // Check that both of the traces made it back to the root r.
  assert(va_trace[na-1] == r);
  assert(vb_trace[nb-1] == r);

  // Count the number of shared edges in the trace.
  int root_length = 0;
  while (va_trace[na-2-root_length] == vb_trace[nb-2-root_length]) {
    ++root_length;
  }

  // If the topo search was performed with efficient terminating conditions,
  // there should be no shared edges.
  assert(root_length == 0);

  // Begin building the cycle.
  *ncycle_out = 0;

  // Add the non-shared nodes of the first trace to the cycle.
  // Include the most recent common ancestor of the two traces.
  for (i=0; i<na-root_length; ++i) {
    cycle_out[(*ncycle_out)++] = va_trace[i];
  }

  // Add the non-shared nodes of the second trace to the cycle.
  // Do not include the most recent common ancestor of the two traces.
  for (i=nb-root_length-2; i>=0; ++i) {
    cycle_out[(*ncycle_out)++] = vb_trace[i];
  }
}


// Get the length of the smallest cycle in the graph,
// and get an associated vertex.
// This vertex may later be used to extract the minimum length cycle.
void get_girth_and_vertex(
    const int *row_ptr, const int *col_ind, int nvertices,
    BFS_WS *bfs_ws, int *depth_ws,
    int *girth_out, int *vertex_out
    )
{
  // If no cycle is found, these values will not change.
  *girth_out = -1;
  *vertex_out = -1;

  // Take the minimum over the minimum cycles observed
  // over all topo sorts over all roots.
  int min_length;
  int root;
  for (root=0; root<nvertices; ++root) {
    int min_length = get_girth_ub(row_ptr, col_ind, nvertices, root,
        bfs_ws, depth_ws);
    if (*girth_out < 0 || min_length < *girth_out) {
      *girth_out = min_length;
      *vertex_out = root;
    }
  }
}


// Get the length of the smallest cycle in the graph.
int get_girth(const int *row_ptr, const int *col_ind, int nvertices,
    BFS_WS *bfs_ws, int *depth_ws
    )
{
  int girth;
  int vertex;
  get_girth_and_vertex(row_ptr, col_ind, nvertices,
      bfs_ws, depth_ws,
      &girth, &vertex);
  return girth;
}


// Get the length of the smallest cycle in the graph,
// and get an associated vertex.
// This vertex may later be used to extract the minimum length cycle.
// This function requires that the graph be connected.
void get_girth_and_vertex_conn(
    const int *row_ptr, const int *col_ind, int nvertices,
    BFS_WS *bfs_ws, int *depth_ws,
    int *girth_out, int *vertex_out
    )
{
  // If no cycle is found, these values will not change.
  *girth_out = -1;
  *vertex_out = -1;

  // Get the number of edges in this graph.
  int nedges = row_ptr[nvertices] - row_ptr[0];

  // If there are not enough edges then the connected graph cannot have a cycle.
  if (nedges < nvertices) {
    return;
  }

  // If there are exactly enough edges then the graph is a pure cycle,
  // so any vertex will work.
  if (nedges == nvertices) {
    *girth_out = nvertices;
    *vertex_out = 0;
    return;
  }

  // Otherwise all cycles of the connected graph
  // will include a vertex of degree greater than two.
  // Look for the minimum such cycle.
  int min_length;
  int root;
  for (root=0; root<nvertices; ++root) {
    int root_degree = row_ptr[root+1] - row_ptr[root];
    if (root_degree > 2) {
      int min_length = get_girth_ub(row_ptr, col_ind, nvertices, root,
          bfs_ws, depth_ws);
      if (*girth_out < 0 || min_length < *girth_out) {
        *girth_out = min_length;
        *vertex_out = root;
      }
    }
  }
}


// Get the length of the smallest cycle in the graph.
// This function requires that the graph be connected.
int get_girth_conn(const int *row_ptr, const int *col_ind, int nvertices,
    BFS_WS *bfs_ws, int *depth_ws
    )
{
  int girth;
  int vertex;
  get_girth_and_vertex_conn(row_ptr, col_ind, nvertices,
      bfs_ws, depth_ws,
      &girth, &vertex);
  return girth;
}

