/*
 * This module is for the grid for the "plan c" solver.
 *
 * The idea is to allocate a 2d grid of integers to track the growing sequence.
 * The first site in the sequence is anchored to the center of the grid,
 * and the size of the allocated grid is large enough so that its
 * radius is greater than the max sequence length.
 *
 * The grid is initially filled with values of -1 indicating lack of presence
 * of a sequence site at the location.
 * As the sequence grows, the value of the grid at sequence site j
 * will have the integer value j.
 * During backtracking, the retracting sequence will refill its formerly
 * occupied grid sites with -1 indicating that the sequence no longer
 * occupies those positions on the grid.
 *
 */

#ifndef PLANCGRID_HEADER
#define PLANCGRID_HEADER

#define RIGHT 0
#define UP 1
#define LEFT 2
#define DOWN 3

#define GRID_EMPTY (-1)
#define GRID_PROBE (-2)
#define GRID_BORDER (-5)


typedef struct tag_GRID
{
  int n;
  int radius;
  int origin_row;
  int origin_col;
  int origin_index;
  int nrows;
  int ncols;
  int area;
  int *data;
} GRID;

void grid_init(GRID *p, int n);

void grid_destroy(GRID *p);

int grid_get(GRID *p, int row, int col);

void grid_set(GRID *p, int row, int col, int value);

void grid_clear(GRID *p, int row, int col);

void print_grid(GRID *p);

#endif

