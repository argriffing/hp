//#include "stdint.h"
//#include "inttypes.h"
//#include "stdlib.h"
//#include "string.h"
//#include "stdbool.h"
//#include "unistd.h"
#include "stdio.h"
#include "assert.h"

#include "plancgrid.h"
#include "compactness.h"


// This is called within the function that evaluates
// a void region for fillability.
void void_init(VOID_INFO *p)
{
  int degree;
  int parity;
  for (parity=0; parity<2; ++parity) {
    for (degree=0; degree<4; ++degree) {
      p->parity_degree_histogram[parity][degree] = 0;
    }
    p->parity_count[parity] = 0;
  }
  p->includes_border = false;
  p->nprobed = 0;
}


// Call this after evaluating the void region for fillability.
void clear_grid_probes(GRID *grid, const int *index_ws, int nprobed)
{
  int i;
  for (i=0; i<nprobed; ++i) {
    int grid_index = index_ws[i];
    int value = grid->data[grid_index];
    if (value != GRID_PROBE) {
      printf("Error: at grid index %d expected %d but found %d\n",
          grid_index, GRID_PROBE, value);
      assert(false);
    }
    grid->data[grid_index] = GRID_EMPTY;
  }
}


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
    )
{
  // Initialize the info about the void region.
  void_init(void_info);

  // Use a breadth first search.
  // Instead of using a csr graph and parent vertices,
  // use the implicit grid graph and the GRID_PROBE markers.
  int i;
  int v;
  int nnext = 1;
  index_ws[0] = void_index;
  assert(grid->data[void_index] == GRID_EMPTY);
  grid->data[void_index] = GRID_PROBE;
  void_info->nprobed++;

  // Define the parity of the last site.
  int parity_of_first_site = (query_index + 1) % 2;
  int parity_of_last_site = (query_index + nremaining) % 2;

  // Define the number of vertices remaining of each parity.
  int nremaining_parity[2];
  nremaining_parity[parity_of_first_site] = (nremaining + 1) / 2;
  nremaining_parity[1 - parity_of_first_site] = nremaining / 2;

  while (nnext) {
    int ncurr = nnext;
    nnext = 0;
    for (i=0; i<ncurr; ++i) {

      // Identify the vertex whose information will be added to the void region
      // info structure and whose neighbors will be added to the
      // next shell of the breadth first search.
      v = index_ws[i];
      assert(v >= 0);
      int parity = (v % 2);

      // Look at the points adjacent to the grid point to be added.
      // The degree of a void space is the number of adjacent void spaces,
      // plus 1 if the query index is adjacent.
      int nborders = 0;
      int degree = 0;
      int direction;
      for (direction=0; direction<4; ++direction) {
        int w = v + delta[direction];
        int value = grid->data[w];
        if (value == GRID_BORDER) {
          printf("found border\n");
          printf("grid index %d\n", v);
          printf("direction %d\n", direction);
          printf("neighbor index %d\n", w);
          printf("neighbor value %d\n", value);
          nborders++;
        } else if (value == GRID_EMPTY) {
          grid->data[w] = GRID_PROBE;
          void_info->nprobed++;
          index_ws[ncurr + nnext] = w;
          nnext++;
          degree++;
        } else if (value == GRID_PROBE || w == query_index) {
          degree++;
        }
      }

      // Require the void degree to be in {1, 2, 3, 4}.
      assert(degree == 1 || degree == 2 || degree == 3 || degree == 4);

      // Require the parity to be in {0, 1}.
      assert(parity == 0 || parity == 1);

      // Add the info corresponding to the current vertex in the void region.
      void_info->parity_count[parity]++;
      void_info->parity_degree_histogram[parity][degree]++;
      if (nborders) {
        void_info->includes_border = true;
      }

      // If the void region includes the border then the region is not fillable.
      if (void_info->includes_border) {
        printf("not fillable: includes border\n");
        int row, col;
        for (row=0; row<grid->nrows; ++row) {
          for (col=0; col<grid->ncols; ++col) {
            printf("%3d ", grid_get(grid, row, col));
          }
          printf("\n");
        }
        return false;
      }

      // If the area of the void region exceeds the number of sites
      // remaining in the sequence then the region is not fillable.
      int pc0 = void_info->parity_count[0];
      int pc1 = void_info->parity_count[1];
      if (pc0 + pc1 > nremaining) {
        printf("not fillable: void area exceeds nremaining\n");
        return false;
      }

      // Among the void spaces with the parity of the last site,
      // at most one space can have degree 1.
      // Among the void spaces with the opposite parity,
      // no space can have degree 1.
      if (void_info->parity_degree_histogram[parity_of_last_site][1] > 1) {
        printf("not fillable: parity of last site degree-1\n");
        return false;
      }
      if (void_info->parity_degree_histogram[1 - parity_of_last_site][1] > 0) {
        printf("not fillable: opposite of parity of last site degree-1\n");
        return false;
      }

      // Does the void contain more spaces of either parity
      // than will be present in the sequence?
      // If so, then the void region cannot be filled.
      int p;
      for (p=0; p<2; ++p) {
        if (void_info->parity_count[p] > nremaining_parity[p]) {
          printf("not fillable: void has too many spaces of parity %d\n", p);
          return false;
        }
      }
    }

    // Skip past the vertices in the previous shell of the flood fill.
    index_ws += ncurr;
  }

  // The region has already been determined to not be adjacent to the border,
  // to not have too many spaces of either parity,
  // and to not have too many void-degree-1 spaces of either parity.
  // Does the void have exactly the right number of spaces of both parities?
  // This will determine whether the region is naively fillable.
  return (
      void_info->parity_count[0] == nremaining_parity[0] &&
      void_info->parity_count[1] == nremaining_parity[1]);
}
