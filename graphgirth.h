// Find the shortest cycle in an undirected unweighted graph.
//
// The length of the shortest cycle is called the girth of the graph.
// The strategy is to do a depth first search from each vertex.
//
// Use the 'compressed row storage' as the graph format.
// http://web.eecs.utk.edu/~dongarra/etemplates/node373.html
// This uses two arrays.
// One array gives the offset corresponding to vertex i.
// The other array gives the sink vertex.
// I will use the somewhat common notation row_ptr and col_ind,
// together with the number of vertices n.
// The degree of vertex v is row_ptr[v+1] - row_ptr[v].


#ifndef GRAPHGIRTH_HEADER
#define GRAPHGIRTH_HEADER

#include "breadthfirst.h"


// Get the length of the smallest cycle detected from topo sort rooted at r.
// This is an upper bound on the girth of the graph.
// The _ws suffixed arrays should be of length nvertices
// and are for temporary storage.
int get_girth_ub(
    const int *row_ptr, const int *col_ind, int r,
    BFS_WS *bfs_ws, int *depth_ws
    );

// Get the smallest cycle detected from topo sort rooted at r.
// The length of this cycle is an upper bound
// on the smallest cycle in the graph.
void get_smallest_cycle_ub(
    const int *row_ptr, const int *col_ind, int r,
    BFS_WS *bfs_ws, int *depth_ws,
    int *va_trace, int *vb_trace,
    int *cycle_out, int *ncycle_out
    );

// Get the length of the smallest cycle in the graph.
int get_girth(const int *row_ptr, const int *col_ind, int nvertices,
    BFS_WS *bfs_ws, int *depth_ws);

// Get the length of the smallest cycle in the graph,
// and get a vertex that is in this cycle.
void get_girth_and_vertex(
    const int *row_ptr, const int *col_ind, int nvertices,
    BFS_WS *bfs_ws, int *depth_ws,
    int *girth_out, int *vertex_out
    );

// Get the length of the smallest cycle in the graph.
// Requires that the graph be connected.
int get_girth_conn(const int *row_ptr, const int *col_ind, int nvertices,
    BFS_WS *bfs_ws, int *depth_ws);

// Get the length of the smallest cycle in the graph,
// and get a vertex that is in this cycle.
// Requires that the graph be connected.
void get_girth_and_vertex_conn(
    const int *row_ptr, const int *col_ind, int nvertices,
    BFS_WS *bfs_ws, int *depth_ws,
    int *girth_out, int *vertex_out
    );

#endif
