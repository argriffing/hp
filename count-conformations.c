#include "stdio.h"
#include "stdint.h"
#include "inttypes.h"
#include "stdlib.h"
#include "string.h"
#include "stdbool.h"
#include "unistd.h"
#include "assert.h"

#include "plancgrid.h"

#define RIGHT 0
#define UP 1
#define LEFT 2
#define DOWN 3

int64_t count_continuations(GRID *grid, int *delta,
    int *direction_histogram, int grid_index, int nsites_remaining)
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
      // With the following restrictions,
      // we compute A046171 instead of A001411.
      if (direction_histogram[RIGHT] == 0 && direction != RIGHT) continue;
      if (direction_histogram[UP] == 0 && direction == DOWN) continue;
      direction_histogram[direction]++;
      int next_grid_index = grid_index + delta[direction];
      if (grid->data[next_grid_index] == GRID_EMPTY) {
        nwalks += count_continuations(grid, delta,
            direction_histogram, next_grid_index, nsites_remaining);
      }
      direction_histogram[direction]--;
    }
  } else {
    nwalks = 1;
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

  // Define a direction histogram which is used for avoiding
  // walks that are redundant due to symmetry.
  int direction_histogram[] = {0, 0, 0, 0};

  // Count the number of walks.
  int64_t nwalks = count_continuations(&grid, delta,
      direction_histogram, grid_index, n);

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
    printf("n %d\n", n);
    printf("walks %" PRId64 "\n", nwalks);
    printf("\n");
    fflush(stdout);
    n++;
  }
  return 0;
}

