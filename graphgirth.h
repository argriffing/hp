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


// Get the length of the smallest cycle detected from topo sort rooted at r.
// This is an upper bound on the girth of the graph.
// The _ws suffixed arrays should be of length nvertices
// and are for temporary storage.
int get_girth_ub(
    const int *row_ptr, const int *col_ind, int nvertices, int r,
    int *parent_ws, int *depth_ws, int *deck_ws, int *next_ws,
    );

// Get the smallest cycle containing vertex r.
// The length of this cycle is an upper bound
// on the smallest cycle in the graph.
int get_smallest_cycle_ub(
    const int *row_ptr, const int *col_ind, int nvertices, int r,
    int *parent_ws, int *depth_ws, int *deck_ws, int *next_ws,
    int *cycle_out, int *ncycle_out
    );

