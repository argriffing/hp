// Breadth first search workspace.
//
// Note that this is a little bit stupid,
// because separate arrays are not needed for deck and next.

#include "stdlib.h"

#include "breadthfirst.h"

void bfs_ws_init(BFS_WS *p, int max_nvertices) {
  int nbytes = max_nvertices * sizeof(int);
  p->parent = (int *) malloc(nbytes);
  p->deck = (int *) malloc(nbytes);
  p->next = (int *) malloc(nbytes);
}

void bfs_ws_destroy(BFS_WS *p) {
  free(p->parent);
  free(p->deck);
  free(p->next);
}

// This is a small breadth first search function.
// The row_ptr and col_ind are parts of a csr graph.
// The parent array should be initialized to -1 at all
// vertices reachable from the set of seed vertices.
// The history array should have a list of seed vertices.
// The nseeds should be the number of seed vertices.
//
// As the breadth first fill progresses,
// the newly visited nodes will be added to the history.
// The return value is the total number of indices
// in the history array, including the initial seed indices.
// 
int bfs_fill(const int *row_ptr, const int *col_ind,
    int *parent, int *history, int nseeds)
{
  int i, j;
  int v, w, seed;
  int ncurr, nnext;
  int ntotal = 0;

  // Initialize parent vertices of seed nodes.
  for (i=0; i<nseeds; ++i) {
    seed = history[i];
    parent[seed] = -2;
  }

  // Do the breadth first fill.
  ncurr = nseeds;
  while (ncurr) {
    nnext = 0;
    for (i=0; i<ncurr; ++i) {
      v = history[i];
      for (j=row_ptr[v]; j<row_ptr[v+1]; ++j) {
        w = col_ind[j];
        if (parent[w] == -1) {
          parent[w] = v;
          history[ncurr + nnext++] = w;
        }
      }
    }
    ntotal += ncurr;
    history += ncurr;
    ncurr = nnext;
  }
  return ntotal;
}

void bfs_clear(int *parent, int *history, int nhistory)
{
  int i, v;
  for (i=0; i<nhistory; ++i) {
    v = history[i];
    parent[v] = -1;
  }
}

