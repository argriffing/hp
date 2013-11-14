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

typedef struct tagS3TABLE {
  int max_ncomponents;
  int max_target_sum;
  int max_table_entries;
  int *data;
} S3TABLE;

void s3table_init(S3TABLE *p, int max_ncomponents, int max_target_sum) {
  p->max_ncomponents = max_ncomponents;
  p->max_target_sum = max_target_sum;
  p->max_table_entries = max_ncomponents * max_target_sum;
  p->data = (int *) malloc(p->max_table_entries * sizeof(int));
}

void s3table_destroy(S3TABLE *p)
{
  free(p->data);
}

void s3table_clear(S3TABLE *p, int ncomponents, int target)
{
  int i, j;
  for (i=0; i<target; ++i) {
    for (j=0; j<ncomponents; ++j) {
      p->data[i*p->max_ncomponents + j] = -1;
    }
  }
}

int s3table_get(S3Table *p, int component, int partial_sum)
{
  return p->data[partial_sum * p->max_ncomponents + component];
}

void s3table_set(S3Table *p, int component, int partial_sum, int contribution)
{
  p->data[partial_sum * p->max_ncomponents + component] = contribution;
}

// The table is required to have been cleared before this call.
// Returns 1 if the target is attainable, otherwise returns 0.
int s3table_forward(S3TABLE *p,
    int *low, int *high, int ncomponents, int target)
{
  int component;
  int partial_sum;
  int contrib;

  // The first column of the table corresponds to the sums
  // that we can reach using only the first component.
  component = 0;
  contrib = 0;
  partial_sum = contrib;
  s3table_set(p, component, partial_sum, contrib);
  for (contrib=low[component]; contrib<=high[component], ++contrib) {
    partial_sum = contrib;
    s3table_set(p, component, partial_sum, contrib);
  }

  // Fill the rest of the columns of the table one by one.
  int prev_comp;
  int prev_part;
  int prev_contrib;
  for (component=1; component<ncomponets-1; ++component)
  {
    int prev_comp = component - 1;

    // When the contribution of the current component is zero,
    // each previously attainable partial sum
    // can be copied to the current column.
    contrib = 0;
    for (prev_part=0; prev_part<=target; ++prev_part) {
      prev_contrib = s3table_get(p, prev_comp, prev_part);
      if (prev_contrib >= 0) {
        s3table_set(p, component, prev_part, contrib)
      }
    }

    // When the contribution of the current component is nonzero,
    // update the attainable partial sums of the current column
    // according to the contribution of the current component.
    for (contrib=low[component]; contrib<=high[component]; ++contrib) {
      for (prev_part=0; prev_part<=target; ++prev_part) {
        prev_contrib = s3table_get(p, prev_comp, prev_part);
        if (prev_contrib >= 0) {
          partial_sum = prev_part + contrib;
          if (partial_sum <= target) {
            s3table_set(p, component, partial_sum, contrib);
          }
        }
      }
    }
  }
  int target_is_attainable = s3table_get(p, ncomponents-1, target) >= 0 ? 1 : 0;
  return target_is_attainable;
}

int s3table_backward(S3TABLE *p,
    int *low, int *high, int ncomponents, int target, int *contribs_out)
{
  // Check that the target is attainable.
  int target_is_attainable = s3table_get(p, ncomponents-1, target) >= 0 ? 1 : 0;
  if (!target_is_attainable) {
    return 0;
  }

  // Trace back through the dynamic programming table,
  // extracting the contribution of each component.
  int current_target = target;
  int component;
  int contrib;
  for (component=ncomponents-1; component>=0; --component) {
    contrib = s3table.get(p, component, target);
    contribs_out[component] = contrib;
    target -= contrib;
  }

  return 1;
}

