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

// Return the new history length, including the original seeds.
int bfs_fill(const int *row_ptr, const int *col_ind,
    int *parent, int *history, int nseeds);

void bfs_clear(int *parent, int *history, int nhistory);

#endif
