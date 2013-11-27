#include "stdio.h"
#include "stdint.h"
#include "inttypes.h"
#include "stdlib.h"
#include "string.h"
#include "stdbool.h"
#include "unistd.h"
#include "assert.h"

#include "plancgrid.h"
#include "compactness.h"

#define RIGHT 0
#define UP 1
#define LEFT 2
#define DOWN 3


///////////////////////////////////////////////////////////////////////////////
// Enumeration of compact self avoiding walks,
// implemented in a way that tries to be a little bit clever.


typedef struct tagSEARCH_OPTION {
  int direction;
  int grid_index;
  bool filling;
} SEARCH_OPTION;

int64_t count_continuations(GRID *grid, const int *delta,
    const int *neighbor_lookup,
    int *direction_histogram, int *index_ws,
    int grid_index, int nsites_remaining, bool filling)
{
  // We are guaranteed that the grid is empty at grid_index.
  // Place the site on the grid and decrement the number of
  // sites remaining to be placed.
  grid->data[grid_index] = 0;
  nsites_remaining--;

  // Explore remaining continuations if possible.
  // If we have placed all of the sites then we are finished
  // with this continuation.
  int64_t nwalks = 0;

  // Check the number of neighboring empty groups.
  int ngroups = count_empty_neighbor_groups(neighbor_lookup,
      grid->data, grid->ncols, grid_index);

  // Declare the search options for continuing the sequence.
  // The extra indirection is for possible sorting of the options,
  // not used here.
  SEARCH_OPTION search_option_data[4];
  SEARCH_OPTION *search_options[4];
  SEARCH_OPTION *p;
  int nsearch_options = 0;
  int i;
  for (i=0; i<4; ++i) {
    search_options[i] = search_option_data + i;
  }

  // Build the array of approved search options.
  int direction;
  for (direction=0; direction<4; ++direction) {
    if (!nsites_remaining) break;
    if (direction_histogram[RIGHT] == 0 && direction != RIGHT) continue;
    if (direction_histogram[UP] == 0 && direction == DOWN) continue;
    if (ngroups >= 3) continue;
    if (filling && ngroups > 1) continue;

    // If the neighboring point is not empty then continue.
    int next_grid_index = grid_index + delta[direction];
    int neighbor = grid->data[next_grid_index];
    if (neighbor != GRID_EMPTY) continue;

    // Check whether the proposed move begins filling a void region.
    // If the void is not fillable then continue.
    bool next_filling = filling;
    if (ngroups == 2 && !filling) {
      if (evaluate_void(grid, delta, index_ws,
          grid_index, next_grid_index, nsites_remaining))
      {
        next_filling = true;
      } else {
        continue;
      }
    }

    // The search option has passed some rigorous screening by this point.
    p = search_options[nsearch_options++];
    p->direction = direction;
    p->grid_index = next_grid_index;
    p->filling = next_filling;
  }

  // Continue exploring the approved search directions if any.
  for (i=0; i<nsearch_options; ++i) {
    p = search_options[i];
    direction_histogram[p->direction]++;
    nwalks += count_continuations(grid, delta,
        neighbor_lookup,
        direction_histogram, index_ws,
        p->grid_index, nsites_remaining, p->filling);
    direction_histogram[p->direction]--;
  }

  // If no sites remain, then count the completed walk
  // unless it introduces a void.
  if (!nsites_remaining) {
    if (ngroups > 1) {
      nwalks = 0;
    } else {
      // debug
      if (!check_compactness(grid, delta, index_ws)) {
        printf("erroneously completed a non-compact walk\n");
        print_grid(grid);
        assert(0);
      }
      nwalks = 1;
    }
  }

  // Roll back the site placement.
  grid->data[grid_index] = GRID_EMPTY;

  // return the number of completed continuations.
  return nwalks;
}


int64_t count_walks(int n)
{
  // Initialize the grid.
  GRID grid;
  grid_init(&grid, n);

  // Define shortcuts for navigating the grid.
  int delta[4];
  delta[LEFT] = -1;
  delta[RIGHT] = 1;
  delta[UP] = -grid.ncols;
  delta[DOWN] = grid.ncols;

  // Define the grid index corresponding to the anchor point.
  int grid_index = grid.origin_row * grid.nrows + grid.origin_col;

  // Index workspace for checking fillability of void regions.
  int index_ws[grid.area];

  // This flag is true if we have begun filling a void region.
  bool filling = false;

  // Initialize neighborhood lookup.
  int neighborhood_lookup[256];
  init_empty_neighbor_group_lookup(neighborhood_lookup);

  int direction_histogram[] = {0, 0, 0, 0};

  // Count the number of walks.
  int64_t nwalks = count_continuations(&grid, delta,
      neighborhood_lookup,
      direction_histogram, index_ws,
      grid_index, n, filling);

  // Destroy the grid.
  grid_destroy(&grid);

  // Return the number of walks.
  return nwalks;
}


///////////////////////////////////////////////////////////////////////////////
// Brute force enumeration of compact self avoiding walks.


int64_t count_continuations_brute(GRID *grid, int *delta,
    int *direction_histogram, int *index_ws,
    int grid_index, int nsites_remaining)
{
  // We are guaranteed that the grid is empty at grid_index.
  // Place the site on the grid and decrement the number of
  // sites remaining to be placed.
  grid->data[grid_index] = 0;
  nsites_remaining--;

  // Fill remaining continuations if possible.
  // If we have placed all of the sites then we are finished
  // with this continuation.
  int64_t nwalks;
  if (nsites_remaining) {
    nwalks = 0;
    int direction;
    for (direction=0; direction<4; ++direction) {
      if (direction_histogram[RIGHT] == 0 && direction != RIGHT) continue;
      if (direction_histogram[UP] == 0 && direction == DOWN) continue;
      direction_histogram[direction]++;
      int next_grid_index = grid_index + delta[direction];
      if (grid->data[next_grid_index] == GRID_EMPTY) {
        nwalks += count_continuations_brute(grid, delta,
            direction_histogram, index_ws,
            next_grid_index, nsites_remaining);
      }
      direction_histogram[direction]--;
    }
  } else {

    // The number of walks is 1 if the conformation is compact,
    // otherwise it is 0.
    nwalks = check_compactness(grid, delta, index_ws);
  }

  // Roll back the site placement.
  grid->data[grid_index] = GRID_EMPTY;

  // return the number of completed continuations.
  return nwalks;
}


// Count the compact walks using brute force.
// Enumerate all walks and throw out the ones that not compact.
int64_t count_walks_brute(int n)
{
  // Initialize the grid.
  GRID grid;
  grid_init(&grid, n);

  // Define shortcuts for navigating the grid.
  int delta[4];
  delta[LEFT] = -1;
  delta[RIGHT] = 1;
  delta[UP] = -grid.ncols;
  delta[DOWN] = grid.ncols;

  // Define the grid index corresponding to the anchor point.
  int grid_index = grid.origin_row * grid.nrows + grid.origin_col;

  // Index workspace for checking compactness;
  int index_ws[grid.area];

  int direction_histogram[] = {0, 0, 0, 0};

  // Count the number of walks.
  int64_t nwalks = count_continuations_brute(&grid, delta, 
      direction_histogram, index_ws, grid_index, n);

  // Destroy the grid.
  grid_destroy(&grid);

  // Return the number of walks.
  return nwalks;
}


///////////////////////////////////////////////////////////////////////////////
// Main function.


int main()
{
  int n=1;
  while (true)
  {
    //int64_t nwalks = count_walks(n);
    int64_t nwalks = count_walks_brute(n);
    printf("n %d\n", n);
    printf("walks %" PRId64 "\n", nwalks);
    printf("\n");
    fflush(stdout);
    n++;
  }
  return 0;
}

