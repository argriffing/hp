#include "assert.h"
#include "stdio.h"

#include "graphgirth.h"


// Get the two nodes furthest away from vertex r
// in the smallest cycle detected in a topo sort rooted at r.
void _girth_ub_helper(
    const int *row_ptr, const int *col_ind, int nvertices, int r,
    int *parent_ws, int *depth_ws, int *deck_ws, int *next_ws,
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
    parent_ws[v] = -1;
    depth_ws[v] = -1;
  }
  parent_ws[r] = -2;
  depth_ws[r] = 0;

  // The root vertex is on deck for breadth first search.
  // The root has depth zero.
  deck_ws[0] = r;
  int ndeck = 1;

  int depth = 0;
  while (ndeck)
  {
    int nnext = 0;
    int ideck;
    for (ideck=0; ideck<ndeck; ++ideck) {
      v = deck_ws[ideck];
      for (i=row_ptr[v]; i<row_ptr[v+1]; ++i) {
        int w = col_ind[i];
        if (w != parent_ws[v]) {

          // Look at each vertex one step outside of the current shell.
          // If the neighbor has already been spotted,
          // then check if the cycle length is the smallest known.
          // If the neighbor has not been spotted
          // then set its parent and depth and add it to the shell.
          if (parent_ws[w] >= 0) {
            int cycle_length = depth_ws[v] + depth_ws[w] + 1;
            if (*min_length < 0 || cycle_length < *min_length) {
              printf("(v, w): (%d, %d)\n", v, w);
              printf("depth of v: %d\n", depth_ws[v]);
              printf("depth of w: %d\n", depth_ws[w]);
              printf("cycle length %d\n", cycle_length);
              printf("\n");
              *min_length = cycle_length;
              *va_out = v;
              *vb_out = w;
            }
          } else {
            parent_ws[w] = v;
            depth_ws[w] = depth_ws[v] + 1;
            next_ws[nnext++] = w;
          }
        }
      }
    }

    // Put the next array on deck.
    int *tmp = deck_ws;
    deck_ws = next_ws;
    next_ws = tmp;
    ndeck = nnext;
  }
}

// Get the length of the smallest cycle detected from topo sort rooted at r.
// This is an upper bound on the girth of the graph.
// The _ws suffixed arrays should be of length nvertices
// and are for temporary storage.
int get_girth_ub(
    const int *row_ptr, const int *col_ind, int nvertices, int r,
    int *parent_ws, int *depth_ws, int *deck_ws, int *next_ws
    )
{
  int min_length, va, vb;
  _girth_ub_helper(row_ptr, col_ind, nvertices, r,
      parent_ws, depth_ws, deck_ws, next_ws,
      &min_length, &va, &vb);
  return min_length;
}

// Get the length of the smallest cycle detected from topo sort rooted at r.
// The length of this cycle is an upper bound
// on the smallest cycle in the graph.
int get_smallest_cycle_ub(
    const int *row_ptr, const int *col_ind, int nvertices, int r,
    int *parent_ws, int *depth_ws, int *deck_ws, int *next_ws,
    int *cycle_out, int *ncycle_out
    )
{
  int min_length, va, vb;
  _girth_ub_helper(row_ptr, col_ind, nvertices, r,
      parent_ws, depth_ws, deck_ws, next_ws,
      &min_length, &va, &vb);

  // Reuse the deck_ws and next_ws buffers for other purposes.
  int *va_trace = deck_ws;
  int *vb_trace = next_ws;

  // Use the parent_ws output from the helper function
  // to trace both vertices back to the root.
  int i, v;
  int na, nb;

  v = va;
  na = 0;
  while (v >= 0) {
    va_trace[na++] = v;
    v = parent_ws[v];
  }

  v = vb;
  nb = 0;
  while (v >= 0) {
    vb_trace[nb++] = v;
    v = parent_ws[v];
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

// Get the length of the smallest cycle and a vertex on the cycle.
int get_girth_and_vertex(
    const int *row_ptr, const int *col_ind, int nvertices
    )
{
  assert(0);
  return 0;
}

// Get the length of the smallest cycle and a vertex on the cycle.
// This function is restricted to connected graphs.
// It treats degree-two vertices more efficiently.
int get_girth_and_vertex_conn(
    const int *row_ptr, const int *col_ind, int nvertices
    )
{
  assert(0);
  return 0;
}

