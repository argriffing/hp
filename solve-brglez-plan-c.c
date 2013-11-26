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
#include "assert.h"

#include "plancgrid.h"
#include "bondsubgraph.h"
#include "sparsetools.h"

#define RIGHT 0
#define UP 1
#define LEFT 2
#define DOWN 3


///////////////////////////////////////////////////////////////////////////////
// These functions are for visualization.


// Fill the conformation string using the current grid.
void fill_conformation_string(char *conformation_string,
    GRID *grid, const int *delta, int n)
{
  // Define the map from directions to letters.
  char compass_rose[] = "xxxx";
  compass_rose[UP] = 'u';
  compass_rose[DOWN] = 'd';
  compass_rose[LEFT] = 'l';
  compass_rose[RIGHT] = 'r';

  // Fill the conformation string.
  int idx = grid->origin_index;
  int v;
  for (v=0; v<n; ++v) {
    assert(idx >= 0);
    int next_idx = -1;
    int direction;
    for (direction=0; direction<4; ++direction) {
      int neighbor_index = idx + delta[direction];
      int neighbor = grid->data[neighbor_index];
      if (neighbor == v+1) {
        next_idx = neighbor_index;
        conformation_string[v] = compass_rose[direction];
      }
    }
    idx = next_idx;
  }
}


// Fill the hp string using the current solution.
void fill_hp_string(char *hp_string, const int *binary_solution, int n)
{
  int v;
  for (v=0; v<n; ++v) {
    char hp_char;
    int value = binary_solution[v];
    if (value == 1) {
      hp_char = 'h';
    } else if (value == 0) {
      hp_char = 'p';
    } else {
      printf("vertex index %d\n", v);
      printf("unrecognized binary solution value %d\n", value);
      assert(0);
    }
    hp_string[v] = hp_char;
  }
}


///////////////////////////////////////////////////////////////////////////////
// Bounds


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
    int nchosen = int_min(degree_histogram[degree], k);
    k -= nchosen;
    max_degree_sum += degree * nchosen;
  }

  // Note that if the conformation is particularly bad,
  // not all of the values of k must be used.
  assert(k >= 0);

  int bound = max_degree_sum / 2;
  return bound;
}


///////////////////////////////////////////////////////////////////////////////
// Solver functions.


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
  GRID *grid;
  int *delta;
  int *distance_from_origin;
  int *degrees;
  int *degree_histogram;
  int *direction_histogram;
  int *binary_solution_ws;
  int *row_ptr;
  int *col_ind;
} PLANC_SOLVER_INFO;



// Compute the potential bond csr graph by tracing the sequence from its origin.
void compute_potential_bond_csr_graph(PLANC_SOLVER_INFO *info)
{
  int i;
  int idx = info->grid->origin_index;
  info->row_ptr[0] = 0;
  for (i=0; i<info->n; ++i) {
    assert(idx > -1);
    info->row_ptr[i+1] = info->row_ptr[i];
    int next_idx = -1;
    int direction;
    for (direction=0; direction<4; ++direction) {
      int neighbor_index = idx + info->delta[direction];
      int neighbor = info->grid->data[neighbor_index];
      if (neighbor != -1) {
        if (neighbor == i+1) {
          next_idx = neighbor_index;
        } else if (abs(i - neighbor) > 1) {
          info->col_ind[info->row_ptr[i+1]++] = neighbor;
        }
      }
    }
    idx = next_idx;
  }

  // Check the csr graph for debugging.
  int csr_fail_flag = 0;
  int nloops = csr_count_loops(info->row_ptr, info->col_ind, info->n);
  if (nloops) {
    printf("found %d loops\n", nloops);
    csr_fail_flag = 1;
  }
  int csr_nvertices = info->row_ptr[info->n] - info->row_ptr[0];
  for (i=0; i<csr_nvertices; ++i) {
    if (info->col_ind[i] < 0 || info->col_ind[i] > info->n - 1) {
      csr_fail_flag = 1;
    }
  }
  if (csr_fail_flag) {
    printf("bad csr graph\n");
    printf("origin index: %d\n", info->grid->origin_index);
    printf("value at origin: %d\n",
        info->grid->data[info->grid->origin_index]);
    print_csr_graph(info->row_ptr, info->col_ind, info->n);
    assert(0);
  }
}


typedef struct tagSEARCH_OPTION {
  int direction;
  int grid_index;
  int distance_from_origin;
} SEARCH_OPTION;

int _compare_search_option_pointers(const void *a, const void *b)
{
  int da = (*(SEARCH_OPTION **)a)->distance_from_origin;
  int db = (*(SEARCH_OPTION **)b)->distance_from_origin;
  return da - db;
}


// Recursive solver.
// TODO remove the recursion.
void rsolve(PLANC_SOLVER_INFO *info,
    int nsteps, int degree_sum, int grid_index)
{
  int v;
  int direction, i;
  int neighbor, neighbor_index, neighbor_degree;

  // We have landed on a grid position that is guaranteed to be empty.
  // Summarize the progress.
  int nsites_curr = nsteps + 1;
  int nsites_remaining = info->n - nsites_curr;

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
  info->degrees[nsteps] = nneighbors;

  // If we are not finished then attempt to continue exploring the sequence.
  if (nsteps < info->n - 1 && info->best_score < info->score_upper_bound) {

    // Construct the search options.
    // Later we will sort the search options so that the points
    // nearer to the origin are explored first.
    SEARCH_OPTION search_option_data[4];
    SEARCH_OPTION *search_options[4];
    SEARCH_OPTION *p;
    int nsearch_options = 0;
    for (i=0; i<4; ++i) {
      search_options[i] = search_option_data + i;
    }

    for (direction=0; direction<4; ++direction) {

      // These conditions avoid exploring redundant conformations.
      if (info->direction_histogram[RIGHT] == 0 && direction != RIGHT) continue;
      if (info->direction_histogram[UP] == 0 && direction == DOWN) continue;

      // If the grid space is empty then add the search option.
      neighbor_index = grid_index + info->delta[direction];
      if (info->grid->data[neighbor_index] == GRID_EMPTY) {
        p = search_options[nsearch_options];
        p->direction = direction;
        p->grid_index = neighbor_index;
        p->distance_from_origin = info->distance_from_origin[p->grid_index];
        nsearch_options++;
      }
    }

    // Sort the search options
    // according to increasing distance from the origin.
    qsort(search_options, nsearch_options,
        sizeof(SEARCH_OPTION *), _compare_search_option_pointers);

    // Explore the search options according to the preferred ordering.
    for (i=0; i<nsearch_options; ++i) {
      p = search_options[i];
      info->direction_histogram[p->direction]++;
      rsolve(info, nsteps+1, degree_sum, p->grid_index);
      info->direction_histogram[p->direction]--;
    }
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
    if (degree_hist_score_bound > info->best_score) {

      // Compute the score using the potential bond csr graph.
      compute_potential_bond_csr_graph(info);
      int score = -1;
      int failflag =  solve_potential_bond_graph(
          info->row_ptr, info->col_ind, info->n, info->k,
          info->binary_solution_ws, &score);
      assert(!failflag);

      // If the score is better than the current best score, then update
      // the best score, the hp string, and the conformation string.
      if (score > info->best_score) {
        info->best_score = score;
        fill_hp_string(info->hp_string,
            info->binary_solution_ws, info->n);
        fill_conformation_string(info->conformation_string,
            info->grid, info->delta, info->n);
      }

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
  info->degrees[nsteps] = -1;
}


// Inputs are sequence length n, sequence weight k, and verbosity.
// Outputs are the conformation string and the hydrophobic/polar string.
int solve(int n, int k, int verbose,
    char *conformation_string, char *hp_string)
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

  // Initialize manhattan distance from origin.
  int distance_from_origin[grid.area];
  int row, col;
  for (row=0; row<grid.nrows; ++row) {
    for (col=0; col<grid.ncols; ++col) {
      int idx = row*grid.ncols + col;
      int d = abs(row - grid.origin_row) + abs(col - grid.origin_col);
      distance_from_origin[idx] = d;
    }
  }

  // Track the potential bond graph degrees
  // of all vertices in the partial sequence.
  int degrees[n];
  memset(degrees, 0, sizeof degrees);

  // Initialize the partial sequence length.
  int nsteps = 0;
  int grid_index = grid.origin_index;
  int degree_histogram[] = {0, 0, 0, 0};
  int direction_histogram[] = {0, 0, 0, 0};
  int degree_sum = 0;

  // Initialize the binary solution workspace for the solver.
  int binary_solution_ws[n];
  memset(binary_solution_ws, -1, sizeof binary_solution_ws);

  // Initialize csr graph workspaces.
  // The row_ptr should have one more space than the number of vertices.
  // The col_ind should have room for twice the max possible
  // number of undirected edges.
  int row_ptr[n+1];
  int col_ind[2*(n+1)];

  // Aggregate solver info.
  // This might look like a lot of things,
  // but it doesn't even include initialization of the 'plan b' solver.
  // The 'plan b' solver is inefficiently initialized and destroyed
  // every time it is used...
  PLANC_SOLVER_INFO info;
  info.n = n;
  info.k = k;
  info.verbose = verbose;
  info.conformation_string = conformation_string;
  info.hp_string = hp_string;
  info.best_score = -1;
  info.score_upper_bound = get_score_upper_bound(n, k);
  info.grid = &grid;
  info.delta = delta;
  info.distance_from_origin = distance_from_origin;
  info.degrees = degrees;
  info.degree_histogram = degree_histogram;
  info.direction_histogram = direction_histogram;
  info.binary_solution_ws = binary_solution_ws;
  info.row_ptr = row_ptr;
  info.col_ind = col_ind;

  // Solve the problem recursively.
  rsolve(&info, nsteps, degree_sum, grid_index);

  grid_destroy(&grid);
  return info.best_score;
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
  int score = solve(n, k, verbose, conformation_string, hp_string);

  // Report the optimal score, as well as the conformation string
  // and hp string for which the optimal score is attained.
  printf("score %d\n", score);
  printf("conformation %s\n", conformation_string);
  printf("hp %s\n", hp_string);

  return 0;
}

