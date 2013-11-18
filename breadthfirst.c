// Breadth first search workspace.

typedef struct BFS_WS {
  int *parent;
  int *deck;
  int *next;
} BFS_WS;

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
