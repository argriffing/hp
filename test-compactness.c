#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "assert.h"

#include "plancgrid.h"
#include "compactness.h"


int _compactness_test_helper(GRID *grid, const int *original_grid,
    int query_index, int void_index, int nremaining,
    bool expected_fillability
    )
{
  int failflag = 0;

  // Define shortcuts for navigating the grid.
  int delta[4];
  delta[LEFT] = -1;
  delta[RIGHT] = 1;
  delta[UP] = -grid->ncols;
  delta[DOWN] = grid->ncols;

  // Stack-allocate the index workspace for flood filling the void regions.
  int index_ws[grid->area];
  memset(index_ws, -1, sizeof index_ws);

  // Evaluate a fillable void space.
  VOID_INFO void_info;
  bool fillable = evaluate_void(grid, delta, &void_info, index_ws,
      query_index, void_index, nremaining);

  // Check the fillability.
  if (fillable != expected_fillability) {
    printf("fail:\n");
    printf("expected fillability %d but observed %d\n",
        expected_fillability, fillable);
    printf("query_index %d\n", query_index);
    printf("void_index %d\n", void_index);
    printf("nremaining %d\n", nremaining);
    printf("\n");
    failflag = 1;
  }

  // Use the grid info to count the number of grid probes.
  int pc0 = void_info.parity_count[0];
  int pc1 = void_info.parity_count[1];
  int nprobed = void_info.nprobed;

  // Clear the grid probes.
  clear_grid_probes(grid, index_ws, nprobed);

  // Check that the grid data has been restored
  // after having cleared the probes.
  int i;
  for (i=0; i<grid->area; ++i) {
    assert(original_grid[i] == grid->data[i]);
  }

  return failflag;
}


int test_compactness_0()
{
  printf("test 0\n");

  // Initialize the test grid.
  int nrows = 7;
  int ncols = 7;
  int area = nrows * ncols;
  int original_grid[] = {
    -5, -5, -5, -5, -5, -5, -5,
    -5, 14, 13, 12, -1, -1, -5,
    -5, 15, -1, 11, -1, -1, -5,
    -5, 16, -1, 10, -1, -1, -5,
    -5, 17, -1, -1, 24, 23, -5,
    -5, 18, 19, 20, 21, 22, -5,
    -5, -5, -5, -5, -5, -5, -5,
  };

  // Initialize the grid structure.
  GRID grid;
  grid_init(&grid, 2);

  // Compare the dimensions to the test grid.
  assert(grid.nrows == nrows);
  assert(grid.ncols == ncols);
  assert(grid.area == area);

  // Copy the test grid data into the grid structure.
  memcpy(grid.data, original_grid, sizeof original_grid);

  // Run some tests.
  int nfails = 0;

  // Declare variables.
  int query_row;
  int query_col;
  int void_row;
  int void_col;
  int query_index;
  int void_index;
  int nremaining;
  int expected_fillability;

  // Setup for vertex 24 stepping to the left.
  query_row = 4;
  query_col = 4;
  void_row = 4;
  void_col = 3;
  query_index = query_row*ncols + query_col;
  void_index = void_row*ncols + void_col;

  // This step should work only when four vertices remain.
  nremaining = 3;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 4;
  expected_fillability = true;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 5;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);

  // Setup for vertex 24 stepping in the up direction.
  query_row = 4;
  query_col = 4;
  void_row = 3;
  void_col = 4;
  query_index = query_row*ncols + query_col;
  void_index = void_row*ncols + void_col;

  // This step should work regardless of how many vertices remain,
  // because this void region is connected to the border.
  nremaining = 3;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 4;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 5;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);

  // Destroy the grid.
  grid_destroy(&grid);

  return nfails;
}


int test_compactness_1()
{
  printf("test 1\n");

  // Initialize the test grid.
  int nrows = 7;
  int ncols = 7;
  int area = nrows * ncols;
  int original_grid[] = {
    -5, -5, -5, -5, -5, -5, -5,
    -5, -1, 10, 10, 10, 10, -5,
    -5, 10, 10, -1, -1, 10, -5,
    -5, 10, -1, -1, -1, 10, -5,
    -5, 10, -1, -1, 10, 10, -5,
    -5, 10, 20, 30, 10, -1, -5,
    -5, -5, -5, -5, -5, -5, -5,
  };

  // Initialize the grid structure.
  GRID grid;
  grid_init(&grid, 2);

  // Compare the dimensions to the test grid.
  assert(grid.nrows == nrows);
  assert(grid.ncols == ncols);
  assert(grid.area == area);

  // Copy the test grid data into the grid structure.
  memcpy(grid.data, original_grid, sizeof original_grid);

  // Run some tests.
  int nfails = 0;

  // Declare variables.
  int query_row;
  int query_col;
  int void_row;
  int void_col;
  int query_index;
  int void_index;
  int nremaining;
  int expected_fillability;

  // Setup for vertex 20 stepping upwards.
  query_row = 5;
  query_col = 2;
  void_row = 4;
  void_col = 2;
  query_index = query_row*ncols + query_col;
  void_index = void_row*ncols + void_col;

  // This step should not work regardless of the number of vertices remaining,
  // because the parity is wrong.
  nremaining = 6;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 7;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 8;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);

  // Setup for vertex 30 stepping in the up direction.
  query_row = 5;
  query_col = 3;
  void_row = 4;
  void_col = 3;
  query_index = query_row*ncols + query_col;
  void_index = void_row*ncols + void_col;

  // This step works when the number of remaining vertices is correct.
  nremaining = 6;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 7;
  expected_fillability = true;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 8;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);

  // Destroy the grid.
  grid_destroy(&grid);

  return nfails;
}


int test_compactness_2()
{
  printf("test 2\n");
  printf("neither starting point is blocked by parity\n");
  printf("but both are blocked by an excess of degree-1 vertices\n");

  // Initialize the test grid.
  int nrows = 7;
  int ncols = 7;
  int area = nrows * ncols;
  int original_grid[] = {
    -5, -5, -5, -5, -5, -5, -5,
    -5, -1, 10, 10, 10, -1, -5,
    -5, 10, 10, -1, 10, -1, -5,
    -5, 10, -1, -1, 10, 10, -5,
    -5, 10, -1, -1, -1, 10, -5,
    -5, 10, 20, 30, 10, 10, -5,
    -5, -5, -5, -5, -5, -5, -5,
  };

  // Initialize the grid structure.
  GRID grid;
  grid_init(&grid, 2);

  // Compare the dimensions to the test grid.
  assert(grid.nrows == nrows);
  assert(grid.ncols == ncols);
  assert(grid.area == area);

  // Copy the test grid data into the grid structure.
  memcpy(grid.data, original_grid, sizeof original_grid);

  // Run some tests.
  int nfails = 0;

  // Declare variables.
  int query_row;
  int query_col;
  int void_row;
  int void_col;
  int query_index;
  int void_index;
  int nremaining;
  int expected_fillability;

  // Setup for vertex 20 stepping upwards.
  query_row = 5;
  query_col = 2;
  void_row = 4;
  void_col = 2;
  query_index = query_row*ncols + query_col;
  void_index = void_row*ncols + void_col;

  // Blocked by too many degree-1 vertices.
  nremaining = 6;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 7;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 8;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);

  // Setup for vertex 30 stepping in the up direction.
  query_row = 5;
  query_col = 3;
  void_row = 4;
  void_col = 3;
  query_index = query_row*ncols + query_col;
  void_index = void_row*ncols + void_col;

  // Blocked by too many degree-1 vertices.
  nremaining = 6;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 7;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);
  //
  nremaining = 8;
  expected_fillability = false;
  nfails += _compactness_test_helper(&grid, original_grid,
      query_index, void_index, nremaining, expected_fillability);

  // Destroy the grid.
  grid_destroy(&grid);

  return nfails;
}

int _help_test_groups(const int *data, int ncols,
    int grid_index, int expected)
{
  printf("help test groups\n");
  int lookup[256];
  init_empty_neighbor_group_lookup(lookup);
  int ngroups = count_empty_neighbor_groups(lookup, data, ncols, grid_index);
  if (ngroups != expected) {
    printf("Error: expected %d empty neighbor groups but saw %d\n",
        expected, ngroups);
    return 1;
  }
  return 0;
}

int test_count_empty_neighbor_groups_0() {
  int data[] = {
    -1, -1, -1,
    -1, -1, -1,
    -1, -1, -1,
  };
  int ncols = 3;
  int grid_index = 4;
  return _help_test_groups(data, ncols, grid_index, 1);
}

int test_count_empty_neighbor_groups_1() {
  int data[] = {
    -1, -1, 10,
    -1, -1, -1,
    -1, 10, -1,
  };
  int ncols = 3;
  int grid_index = 4;
  return _help_test_groups(data, ncols, grid_index, 2);
}

int test_count_empty_neighbor_groups_2() {
  int data[] = {
    10, -1, -1,
    -1, -1, -1,
    -1, -1, -1,
  };
  int ncols = 3;
  int grid_index = 4;
  return _help_test_groups(data, ncols, grid_index, 1);
}

int test_count_empty_neighbor_groups_3() {
  int data[] = {
    10, -1, 10,
    -1, 42, -1,
    10, -1, 10,
  };
  int ncols = 3;
  int grid_index = 4;
  return _help_test_groups(data, ncols, grid_index, 4);
}

int test_count_empty_neighbor_groups_4() {
  int data[] = {
    -1, 10, -1,
    -1, 42, 10,
    -1, -1, -1,
  };
  int ncols = 3;
  int grid_index = 4;
  return _help_test_groups(data, ncols, grid_index, 1);
}


int main()
{
  int nfails = 0;

  printf("testing compact neighbor group counts...\n");
  nfails += test_count_empty_neighbor_groups_0();
  nfails += test_count_empty_neighbor_groups_1();
  nfails += test_count_empty_neighbor_groups_2();
  nfails += test_count_empty_neighbor_groups_3();
  nfails += test_count_empty_neighbor_groups_4();

  printf("testing fillability...\n");
  nfails += test_compactness_0();
  nfails += test_compactness_1();
  nfails += test_compactness_2();

  if (nfails) {
    printf("failed testing: %d tests failed\n", nfails);
  } else {
    printf("success: all tests passed\n");
  }
}

