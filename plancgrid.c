#include "plancgrid.h"

#include "stdlib.h"


void grid_init(GRID *p, int n)
{
  p->n = n;
  p->radius = n;
  p->nrows = 2*p->radius + 1;
  p->ncols = 2*p->radius + 1;
  p->area = p->nrows * p->ncols;
  p->origin_row = p->radius;
  p->origin_col = p->radius;
  p->origin_index = p->origin_row * p->ncols + p->origin_col;
  p->data = (int *) malloc(p->area * sizeof(int));
  int i;
  for (i=0; i<p->area; ++i) {
    p->data[i] = -1;
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

