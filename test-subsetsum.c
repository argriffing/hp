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
  s3table_clear(&t, ncomponents, target);
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
  s3table_destroy(&t);

  return failflag;
}

int t_a0()
{
  char name[] = "first test";
  int low[] = {10, 6, 12};
  int high[] = {10, 8, 12};
  int expected[] = {10, 0, 12};
  int ncomponents = 3;
  int target = 22;
  return test(name, low, high, expected, ncomponents, target);
}

int t_a1()
{
  char name[] = "check first of three components";
  int low[] = {10, 6, 12};
  int high[] = {10, 8, 12};
  int expected[] = {10, 0, 0};
  int ncomponents = 3;
  int target = 10;
  return test(name, low, high, expected, ncomponents, target);
}

int t_a2()
{
  char name[] = "check second of three components";
  int low[] = {10, 6, 12};
  int high[] = {10, 8, 12};
  int expected[] = {0, 6, 0};
  int ncomponents = 3;
  int target = 6;
  return test(name, low, high, expected, ncomponents, target);
}

int t_a3()
{
  char name[] = "check third of three components";
  int low[] = {10, 6, 12};
  int high[] = {10, 8, 12};
  int expected[] = {0, 0, 12};
  int ncomponents = 3;
  int target = 12;
  return test(name, low, high, expected, ncomponents, target);
}

int t1()
{
  char name[] = "single component sum zero attainable";
  int low[] = {3};
  int high[] = {6};
  int expected[] = {0};
  int ncomponents = 1;
  int target = 0;
  return test(name, low, high, expected, ncomponents, target);
}

int t2()
{
  char name[] = "single component sum two unattainable";
  int low[] = {3};
  int high[] = {6};
  int expected[] = {-1};
  int ncomponents = 1;
  int target = 2;
  return test(name, low, high, expected, ncomponents, target);
}

int t3()
{
  char name[] = "single component sum five attainable";
  int low[] = {3};
  int high[] = {6};
  int expected[] = {5};
  int ncomponents = 1;
  int target = 5;
  return test(name, low, high, expected, ncomponents, target);
}

int main()
{
  int nfails = 0;

  nfails += t_a0();
  nfails += t_a1();
  nfails += t_a2();
  nfails += t_a3();
  nfails += t1();
  nfails += t2();
  nfails += t3();

  if (nfails) {
    printf("failed testing: %d tests failed\n", nfails);
  } else {
    printf("success: all tests passed\n");
  }
}

