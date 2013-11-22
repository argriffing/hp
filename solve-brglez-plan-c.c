// Build with gcc using a command like
// $ gcc -O3 -march=native solve-brglez-plan-c.c -o solve-brglez-plan-c
//
// Call like
// ./solve-brglez-plan-c -n 22 -k 10
//
// Profile using these commands.
// $ valgrind --tool=callgrind ./a.out <args>
// $ kcachegrind callgrind.out.x

#include "stdio.h"
#include "stdint.h"
#include "stdlib.h"
#include "string.h"
#include "stdbool.h"
#include "unistd.h"

#include "plancgrid.h"

#define RIGHT 0
#define UP 1
#define LEFT 2
#define DOWN 3

int int_max(int a, int b) {return a > b ? a : b;}
int int_min(int a, int b) {return a < b ? a : b;}


// These bounds are hardcoded according to mathematical analysis,
// for example using the "parity" upper bound.
int get_score_upper_bound(int n, int k)
{
  if (k < 2) {
    return 0;
  } else if (k == 2) {
    return 1;
  } else if (k == 3) {
    return 2;
  } else if (k == 4) {
    return 4;
  } else if (k == 5) {
    return 5;
  } else if (k == 7) {
    return 7;
  } else if ((n - k) % 2) {
    return k;
  } else {
    return k+1;
  }
}


// Given the current degree sum and the number of sites remaining,
// get an upper bound on the score.
// This bound applies to partial sequences.
int get_degree_sum_score_bound(int degree_sum, int nsites_remaining)
{
  int max_degree_sum = degree_sum;
  if (nsites_remaining) {
    max_degree_sum += 2 * (2*nsites_remaining + 1);
  }
  int bound = max_degree_sum / 2;
  return bound;
}


// Given the numbers of vertices with degrees 0, 1, 2, and 3,
// get the score upper bound on the number of bonds that could be formed
// by taking a subgraph of k vertices.
// This bound applies only to the full sequence, not to partial sequences.
// TODO this could be improved by distinguishing between vertex parities.
int get_degree_hist_score_bound(const int *degree_histogram, int k)
{
  int max_degree_sum = 0;
  int degree;
  for (degree=3; degree>=0; --degree) {
    int nchosen = int_min(histogram[degree], k);
    k -= nchosen;
    max_degree_sum += degree * nchosen;
  }
  assert(k == 0);
  int bound = max_degree_sum / 2;
  return bound;
}


// This structure just aggregates some things
// that do not change within the recursive solver.
// It is not really a coherent object.
typedef struct tagPLANC_SOLVER_INFO {
  int n;
  int k;
  int verbose;
  char *conformation_string;
  char *hp_string;
  int best_score;
  int score_upper_bound;
  GRID *pgrid;
  int *delta;
  int *degrees;
  int *degree_histogram;
} PLANC_SOLVER_INFO;


// Recursive solver.
// TODO remove the recursion.
void rsolve(PLANC_SEARCH_INFO *info,
    int nsteps, int degree_sum, int grid_index)
{
  int direction, i;
  int neighbor, neighbor_index, neighbor_degree;

  // We have landed on a grid position that is guaranteed to be empty.
  // Summarize the progress.
  int nsites_curr = nsteps + 1;
  int nsites_remaining = n - nsites_curr;

  // Make a list of neighbors that are not adjacent along the primary sequence.
  int neighbors[3];
  int nneighbors = 0;
  for (direction=0; direction<4; ++direction) {
    neighbor_index = grid_index + info->delta[direction];
    neighbor = info->grid->data[neighbor_index];
    if (neighbor != -1 && abs(nsteps - neighbor) > 1) {
      neighbors[nneighbors++] = neighbor;
    }
  }

  // Update the degree sum so that we can check a score bound.
  degree_sum += 2 * nneighbors;

  // Check the score bound that uses the degree sum.
  // If we have no hope of beating the best score
  // then do not explore the continuations of this partial sequence.
  int degree_sum_score_bound = get_degree_sum_score_bound(
      degree_sum, nsites_remaining);
  if (degree_sum_score_bound <= info->best_score) {
    return;
  }

  // Update some quantities that will need to be rolled back
  // before we return from this function. This includes:
  // - update the grid
  // - update the array of vertex degrees
  // - update the degree histogram
  // These could theoretically be handled by the stack,
  // but we use mutability because it is more efficient.
  info->grid->data[grid_index] = nsteps;
  for (i=0; i<nneighbors; ++i) {
    neighbor = neighbors[i];
    neighbor_degree = info->degrees[neighbor];
    info->degrees[neighbor]++;
    info->degree_histogram[neighbor_degree]--;
    info->degree_histogram[neighbor_degree+1]++;
  }
  info->degree_histogram[nneighbors]++;

  // Check the sequence continuation in each direction.
  // For the first two steps, the number of directions is limited
  // to avoid redundantly searching conformations
  // that are equivalent by symmetry.
  direction = 0;
  while (info->best_score < info->score_upper_bound) {
    if (nsteps == info->n - 1) break;
    if (nsteps == 0 && direction > 0) break;
    if (nsteps == 1 && direction > 1) break;
    if (direction > 3) break;

    // Compute the grid index of the next site in the sequence.
    // If the grid is not empty at that position,
    // then continue to the next direction.
    neighbor_index = grid_index + info->delta[direction];
    if (info->grid->data[neighbor_index] != -1) {
      continue;
    }

    // The grid is empty in the direction of interest,
    // so proceed with the recursion.
    rsolve(info, nsteps+1, degree_sum, neighbor_index);

    ++direction;
  }

  // If the sequence is complete then compute a tighter score upper bound
  // using the degree histogram.
  // If the upper bound is better than the best score,
  // then fully evaluate the score of the complete sequence.
  // If the fully evaluated score is better than the best score,
  // then update the search summary.
  if (nsteps == info->n - 1) {

    // Compute an upper bound.
    int degree_hist_score_bound = get_degree_hist_score_bound(
        info->degree_histogram, info->k);
    if (degree_hist_score_bound > best_score) {

      // Compute the potential bond csr graph
      // by tracing the sequence from its origin.
      int idx = info->grid->origin_index;
      int next_idx;
      info->row_ptr[0] = 0;
      for (i=0; i<info->n; ++i) {
        next_idx = -1;
        for (direction=0; direction<4; ++direction) {
          neighbor_index = idx + info->delta[direction];
          neighbor = info->grid->data[neighbor_index];
          if (neighbor != -1) {
            if (neighbor == i+1) {
              next_idx = neighbor_index;
            } else if (abs(i - neighbor) > 1) {
              info->col_ind[row_ptr[i+1]++] = neighbor;
            }
          }
        }
        idx = next_idx;
      }

      // Compute the score using the potential bond csr graph.

    } // end degree hist upper bound check
  }

  // Roll back.
  info->grid->data[grid_index] = -1;
  for (i=0; i<nneighbors; ++i) {
    neighbor = neighbors[i];
    neighbor_degree = info->degrees[neighbor];
    info->degrees[neighbor]--;
    info->degree_histogram[neighbor_degree]--;
    info->degree_histogram[neighbor_degree-1]++;
  }
  info->degree_histogram[nneighbors]--;
}


// Inputs are sequence length n, sequence weight k, and verbosity.
// Outputs are the conformation string and the hydrophobic/polar string.
void solve(int n, int k, int verbose,
    char *conformation_string, char *hp_string)
{
  int i;
  int degree_histogram[] = {0, 0, 0, 0};
  int degree_sum = 0;

  // Initialize the grid.
  GRID grid;
  grid_init(&grid, n);

  // Define shortcuts for navigating the grid.
  int delta[4];
  delta[LEFT] = -1;
  delta[RIGHT] = 1;
  delta[UP] = -grid.ncols;
  delta[DOWN] = grid.ncols;

  // Track the potential bond graph degrees
  // of all vertices in the partial sequence.
  int degrees[n];
  memset(degrees, 0, sizeof degrees);

  // Track the best score.
  // Negative score means that no feasible solution has been found.
  int best_score = -1;

  // Initialize the partial sequence length.
  int nsteps = 0;

  // Put the initial vertex at the origin of the grid.
  current_grid_index = grid.origin_row * grid.nrows + grid.origin_col;
  grid.data[current_grid_index] = 0;

  // Do the search.
  // Stop the search if we reach the theoretical upper bound.
  // Also stop the search if we exhaust the search space.
  int score_upper_bound = get_score_upper_bound(n, k);
  while (0 <= nsteps && best_score < score_upper_bound) {

    int ncurr = nsteps + 1;
    int nsites_remaining = n - ncurr;

    // Check the score bound that uses the degree sum.
    // If we have no hope of beating the best score then backtrack.
    int degree_sum_score_bound = get_degree_sum_score_bound(
        degree_sum, nsites_remaining);
    if (degree_sum_score_bound <= best_score) {
      ncurr--;
      continue;
    }

    if (!nsites_remaining) {

      // Check the score bound that uses the degree histogram.
      // If we have no hope of beating the best score then backtrack.
      int degree_hist_score_bound = get_degree_hist_score_bound(
          degree_histogram, k);
      if (degree_sum_score_bound <= best_score) {
        ncurr--;
        continue;
      }

      // Fully evaluate the complete sequence.
      ;
    }

    // At this point we have an incomplete sequence
    // whose completion will not obviously fail to beat the best score.
    // If the direction is positive then we have already failed at
    // attempting a search direction
    // 4 then we have tried all of the directions.
    // Otherwise try the current direction
  }

  grid_destroy(&grid);
  return best_score;
}


int main(int argc, char **argv)
{

  // Read sequence length n and sequence weight k from the command line.
  // Also read the verbosity option.
  int c;
  int n = -1;
  int k = -1;
  char *vopt = 0;
  char *nopt = 0;
  bool verbose = false;
  while ( (c = getopt(argc, argv, "vn:k:")) != -1)
  {
    int this_option_optind = optind ? optind : 1;
    switch (c)
    {
      case 'v':
        verbose = true;
        break;
      case 'n':
        n = atoi(optarg);
        break;
      case 'k':
        k = atoi(optarg);
        break;
      default:
        printf("?? getopt returned character code 0%o ??\n", c);
    }
  }

  if (n < 4)
  {
    printf("Error: sequence length n must be at least 4\n");
    return 1;
  }

  if (k < 2)
  {
    printf("Error: sequence weight k must be at least 2\n");
    return 1;
  }

  // Report the inputs.
  printf("n %d\n", n);
  printf("k %d\n", k);

  // Allocate room for the conformation string.
  // This needs n-1 letters plus the null termination.
  // The letters are in {l, u, r, d} for 'left', 'up', 'right', 'down'.
  char conformation_string[n];
  memset(conformation_string, 0, sizeof conformation_string);

  // Allocate room for the hydrophobi/polar string.
  // This needs n letters plus the null termination.
  // The letters are in {h, p} for 'hydrophobic', 'polar'.
  char hp_string[n+1];
  memset(hp_string, 0, sizeof hp_string);

  // Solve the puzzle.
  int score = solve(n, k, v, conformation_string, hp_string);

  // Report the optimal score, as well as the conformation string
  // and hp string for which the optimal score is attained.
  printf("score %d\n", score);
  printf("conformation %s\n", conformation_string);
  printf("hp %s\n", hp_string);

  return 0;
}

