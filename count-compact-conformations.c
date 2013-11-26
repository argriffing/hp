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
  int direction;
  int next_grid_index;

  // Check the number of neighboring empty groups.
  int ngroups = count_empty_neighbor_groups(neighbor_lookup,
      grid->data, grid->ncols, grid_index);

  if (nsites_remaining) {

    // If the number of neighboring empty groups is less than two
    // then carry on as usual.
    if (ngroups < 2) {
      for (direction=0; direction<4; ++direction) {
        if (direction_histogram[RIGHT] == 0 && direction != RIGHT) continue;
        if (direction_histogram[UP] == 0 && direction == DOWN) continue;
        direction_histogram[direction]++;
        next_grid_index = grid_index + delta[direction];
        if (grid->data[next_grid_index] == GRID_EMPTY) {
          nwalks += count_continuations(grid, delta,
              neighbor_lookup,
              direction_histogram, index_ws,
              next_grid_index, nsites_remaining, filling);
        }
        direction_histogram[direction]--;
      }
    }
    
    // If the number of groups is equal to two
    // and we are not current filling,
    // then explore only the fillable directions if any,
    // setting the filling state for the remaining continuations.
    if (ngroups == 2 && !filling) {
      for (direction=0; direction<4; ++direction) {
        if (direction_histogram[RIGHT] == 0 && direction != RIGHT) continue;
        if (direction_histogram[UP] == 0 && direction == DOWN) continue;
        direction_histogram[direction]++;
        next_grid_index = grid_index + delta[direction];
        if (grid->data[next_grid_index] == GRID_EMPTY) {
          bool next_fillable = evaluate_void(
              grid, delta, index_ws,
              grid_index, next_grid_index, nsites_remaining);
          if (next_fillable) {
            nwalks += count_continuations(grid, delta,
                neighbor_lookup,
                direction_histogram, index_ws,
                next_grid_index, nsites_remaining, next_fillable);
          }
        }
        direction_histogram[direction]--;
      }
    }
    
    // If the number of groups is equal to two and we are filling,
    // or if the number of groups is greater than two,
    // then we do not explore that direction.
  } else {
    if (ngroups > 1) {
      nwalks = 0;
    } else {
      // debug
      //if (!check_compactness(grid, delta, index_ws)) {
        //printf("erroneously completed a non-compact walk\n");
        //print_grid(grid);
        //assert(0);
      //}
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


int main()
{
  int n=1;
  while (true)
  {
    int64_t nwalks = count_walks(n);
    //int64_t nwalks = count_walks_brute(n);
    printf("n %d\n", n);
    printf("walks %" PRId64 "\n", nwalks);
    printf("\n");
    fflush(stdout);
    n++;
  }
  return 0;
}

