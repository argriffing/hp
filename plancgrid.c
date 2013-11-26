#include "plancgrid.h"

#include "stdio.h"
#include "stdlib.h"


void grid_init(GRID *p, int n)
{
  p->n = n;
  p->radius = n+1;
  p->nrows = 2*p->radius + 1;
  p->ncols = 2*p->radius + 1;
  p->area = p->nrows * p->ncols;
  p->origin_row = p->radius;
  p->origin_col = p->radius;
  p->origin_index = p->origin_row * p->ncols + p->origin_col;
  p->data = (int *) malloc(p->area * sizeof(int));

  // Clear the grid.
  int i;
  for (i=0; i<p->area; ++i) {
    p->data[i] = GRID_EMPTY;
  }

  // Add the borders.
  int row, col;
  for (row=0; row<p->nrows; ++row) {
    grid_set(p, row, 0, GRID_BORDER);
    grid_set(p, row, p->ncols-1, GRID_BORDER);
  }
  for (col=0; col<p->ncols; ++col) {
    grid_set(p, 0, col, GRID_BORDER);
    grid_set(p, p->nrows-1, col, GRID_BORDER);
  }
}

void print_grid(GRID *p)
{
  int row, col;
  printf("%d rows  %d cols\n", p->nrows, p->ncols);
  for (row=0; row<p->nrows; ++row) {
    for (col=0; col<p->ncols; ++col) {
      printf("%3d ", grid_get(p, row, col));
    }
    printf("\n");
  }
}

void grid_destroy(GRID *p)
{
  free(p->data);
}

int grid_get(GRID *p, int row, int col)
{
  return p->data[row * p->ncols + col];
}

void grid_set(GRID *p, int row, int col, int value)
{
  p->data[row * p->ncols + col] = value;
}

void grid_clear(GRID *p, int row, int col)
{
  grid_set(p, row, col, -1);
}

