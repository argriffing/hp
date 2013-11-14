#include "subsetsum.h"

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

int s3table_get(S3TABLE *p, int component, int partial_sum)
{
  return p->data[partial_sum * p->max_ncomponents + component];
}

int s3table_attainable(S3TABLE *p, int component, int partial_sum)
{
  return s3table_get(p, component, partial_sum) >= 0;
}

void s3table_set(S3TABLE *p, int component, int partial_sum, int contribution)
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
      if (s3table_attainable(p, prev_comp, prev_part)) {
        s3table_set(p, component, prev_part, contrib)
      }
    }

    // When the contribution of the current component is nonzero,
    // update the attainable partial sums of the current column
    // according to the contribution of the current component.
    for (contrib=low[component]; contrib<=high[component]; ++contrib) {
      for (prev_part=0; prev_part<=target; ++prev_part) {
        if (s3table_attainable(p, prev_comp, prev_part)) {
          partial_sum = prev_part + contrib;
          if (partial_sum <= target) {
            s3table_set(p, component, partial_sum, contrib);
          }
        }
      }
    }
  }

  return s3table_attainable(p, ncomponents-1, target);
}

// The table is expected to have been processed by a forward pass.
// Returns 1 if the target is attainable, otherwise returns 0.
int s3table_backward(S3TABLE *p,
    int ncomponents, int target, int *contribs_out)
{
  int component;

  // If the target is not attainable then return 0
  // and fill the contribs with -1.
  if (!s3table_attainable(p, ncomponents-1, target)) {
    for (component=0; component<ncomponents; ++component) {
      contribs_out[component] = -1;
    }
    return 0;
  }

  // Trace back through the dynamic programming table,
  // extracting the contribution of each component.
  int current_target = target;
  int contrib;
  for (component=ncomponents-1; component>=0; --component) {
    contrib = s3table.get(p, component, target);
    contribs_out[component] = contrib;
    target -= contrib;
  }

  return 1;
}

