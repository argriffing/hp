#include "stdio.h"
#include "stdlib.h"

#include "subsetsum.h"

#define PASS 0
#define FAIL 1

// Return 1 if failed.
int test(const char *test_name,
    int *low, int *high, int *expected, int ncomponents, int target)
{
  int *contribs = (int *) calloc(ncomponents, sizeof(int));
  S3TABLE t;
  s3table_init(&t, ncomponents, target);
  s3table_forward(&t, low, high, ncomponents, target);
  s3table_backward(&t, ncomponents, target, contribs);
  int failflag = 0;
  int i;
  for (i=0; i<ncomponents; ++i) {
    if (expected[i] != contribs[i]) {
      failflag = 1;
    }
  }

  printf("failed a test\n");
  printf("test name: %s\n", test_name);
  
  printf("target sum: %d\n", target);

  printf("low component values: ");
  for (i=0; i<ncomponents; ++i) {
    printf("%4d ", low[i]);
  }

  printf("high component values: ");
  for (i=0; i<ncomponents; ++i) {
    printf("%4d ", high[i]);
  }

  printf("expected component contribs: ");
  for (i=0; i<ncomponents; ++i) {
    printf("%4d ", expected[i]);
  }

  printf("observed component contribs: ");
  for (i=0; i<ncomponents; ++i) {
    printf("%4d ", expected[i]);
  }

  free(contribs);
  return failflag;
}

int t0()
{
  int low[] = {10, 6, 12};
  int high[] = {10, 8, 12};
  int expected[] = {1, 0, 1}
  int ncomponents = 3;
  int target = 21;
  return test(low, high, expected, ncomponets, target);
}

int main()
{
  int nfails = 0;

  nfails += t0();

  if (nfails) {
    printf("failed testing: %d tests failed", nfails);
  } else {
    printf("success: all tests passed\n");
  }
}

