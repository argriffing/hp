// Breadth first search workspace.

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

