```
$ gcc -o solve-brglez-plan-b solve-brglez-plan-b.c sparsetools.c subsetsum.c \
graphgirth.c connectedcomponents.c breadthfirst.c bondsubgraph.c
$ echo dlluruluurdrdrurddl | ./solve-brglez-plan-b 10
9
10100110100101100101
```

Here's an example of constructing and solving a "plan b" 2d hp puzzle
from the example subset sum problem instance from
[here](http://www.math.sunysb.edu/~scott/blair/Subset_sum_problems_are.html).
The target sum in the example is 5842 which is multiplied by 4
because of technical details in the problem reduction.
The transformed problem has over 70,000 vertices,
and by finding a coloring that has one bond per sequence weight,
it solves the unary subset sum problem whose target sum
is 1/4 the sum of the requested sequence weight.

```
$ python hp-subset-sum.py 267 493 869 961 1000 1153 1246 1598 1766 1922 > puzzle.out
$ ./solve-brglez-plan-b 23368 < puzzle.out | head -n 1
23368
```
