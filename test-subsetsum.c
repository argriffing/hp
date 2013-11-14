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

  if (failflag)
  {

    printf("failed a test\n");
    printf("test name: %s\n", test_name);
    
    printf("target sum: %d\n", target);

    printf("component values: ");
    for (i=0; i<ncomponents; ++i) {
      if (low[i] == high[i]) {
        printf("%d ", low[i]);
      } else {
        printf("%d..%d ", low[i], high[i]);
      }
    }
    printf("\n");

    printf("expected component contribs: ");
    for (i=0; i<ncomponents; ++i) {
      printf("%4d ", expected[i]);
    }
    printf("\n");

    printf("observed component contribs: ");
    for (i=0; i<ncomponents; ++i) {
      printf("%4d ", contribs[i]);
    }
    printf("\n");

    printf("\n");
  }

  free(contribs);
  return failflag;
}

int t0()
{
  char name[] = "first test";
  int low[] = {10, 6, 12};
  int high[] = {10, 8, 12};
  int expected[] = {1, 0, 1};
  int ncomponents = 3;
  int target = 22;
  return test(name, low, high, expected, ncomponents, target);
}

int main()
{
  int nfails = 0;

  nfails += t0();

  if (nfails) {
    printf("failed testing: %d tests failed\n", nfails);
  } else {
    printf("success: all tests passed\n");
  }
}

