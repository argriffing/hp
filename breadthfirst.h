// Breadth first search workspace.

#ifndef BREADTHFIRST_HEADER
#define BREADTHFIRST_HEADER

typedef struct BFS_WS {
  int *parent;
  int *deck;
  int *next;
} BFS_WS;

void bfs_ws_init(BFS_WS *p, int max_nvertices);

void bfs_ws_destroy(BFS_WS *p);

#endif
