#ifndef COMPACTNESS_HEADER
#define COMPACTNESS_HEADER

#include "stdbool.h"

#include "plancgrid.h"




int count_empty_neighbor_groups(const int *lookup,
    const int *data, int ncols, int grid_index);

void init_empty_neighbor_group_lookup(int *lookup);

// Detect properties of a connected region of empty grid points.
// The query_index is the grid index of the sequence site
// which is adjacent to the void.
// The void_index is the grid index of the void site contained
// in the void region.
// When the width of the grid is odd, we can just use the
// parity of the grid index.
//
// Return true if the void is fillable, otherwise return false.
//
bool evaluate_void(GRID *grid, const int *delta, int *index_ws,
    int query_index, int void_index, int nremaining);

int check_compactness(GRID *grid, const int *delta, int *index_ws);

#endif

