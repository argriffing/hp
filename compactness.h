#ifndef COMPACTNESS_HEADER
#define COMPACTNESS_HEADER

#include "stdbool.h"

#include "plancgrid.h"



typedef struct tagVOID_INFO {
  int parity_degree_histogram[2][5];
  int parity_count[2];
  bool includes_border;
  int nprobed;
} VOID_INFO;


int count_empty_neighbor_groups(const int *lookup,
    const int *data, int ncols, int grid_index);


// This is called within the function that evaluates
// a void region for fillability.
void void_init(VOID_INFO *p);


// Call this after evaluating the void region for fillability.
void clear_grid_probes(GRID *grid, const int *index_ws, int nprobed);


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
bool evaluate_void(GRID *grid, const int *delta,
    VOID_INFO *void_info, int *index_ws,
    int query_index, int void_index, int nremaining
    );

#endif

