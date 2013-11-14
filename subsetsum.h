// Subset sum.
//
// This solves a weird version of the subset sum problem,
// using the pseudo-polynomial dynamic programming solution.
// In particular, is implements the dynamic programming algorithm
// referred to in the last sentence of step 5 of the proof of proposition 4
// in "Complexity of finding dense subgraphs" by Asahiro et al. 2002.
//
// In this variant, the input is a set of components and a target sum.
// Each component c_i is allowed to contribute some amount to the sum.
// A component can contribute either 0 (if not chosen)
// or it can contribute some amount between c_i_low and c_i_high inclusive.
// The puzzle is to pick which components to use,
// and how much each component c_i contributes,
// to exactly attain the target sum if possible.

// Include guards.
#ifndef SUBSET_HEADER
#define SUBSET_HEADER

typedef struct tagS3TABLE {
  int max_ncomponents;
  int max_target_sum;
  int max_table_entries;
  int *data;
} S3TABLE;

void s3table_init(S3TABLE *p, int max_ncomponents, int max_target_sum);

void s3table_destroy(S3TABLE *p);

void s3table_clear(S3TABLE *p, int ncomponents, int target);

int s3table_get(S3TABLE *p, int component, int partial_sum);

void s3table_attainable(S3TABLE *p, int component, int partial_sum);

void s3table_set(S3TABLE *p, int component, int partial_sum, int contribution);

// The table is required to have been cleared before this call.
// Returns 1 if the target is attainable, otherwise returns 0.
int s3table_forward(S3TABLE *p,
    int *low, int *high, int ncomponents, int target);

// The table is expected to have been processed by a forward pass.
// Returns 1 if the target is attainable, otherwise returns 0.
int s3table_backward(S3TABLE *p,
    int *low, int *high, int ncomponents, int target, int *contribs_out);

#endif

