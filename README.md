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
because of technical details in the problem reduction.  The transformed
[problem](https://raw.github.com/argriffing/hp/master/puzzle.out)
has over 70,000 vertices, and by finding a coloring
that has one bond per sequence weight,
it solves the unary subset sum problem whose target sum
is 1/4 the sum of the requested sequence weight.

```
$ python hp-subset-sum.py 267 493 869 961 1000 1153 1246 1598 1766 1922 > puzzle.out
$ ./solve-brglez-plan-b 23368 < puzzle.out | head -n 1
23368
```

oeis sequences
--------------

The sequence of counts of rooted self avoiding walks is
[A001411](https://oeis.org/A001411)
and related sequences can be found by searching for
[references](https://oeis.org/search?q=A001411)
to this sequence in the database.
For example the sequence
[A046171](https://oeis.org/A046171)
counts rooted self avoiding walks in a way that removes
much of the redundancy caused by symmetry.


notes on various definitions of compactness
-------------------------------------------

In
Predicting the non-compact conformation of amino acid sequence by
particle swarm optimization
by Yuzhen Guo and Yong Wang (2013)
doi 10.1109/ISB.2013.6623805
compactness is defined as follows.
Suppose that the number of amino acids in a protein sequence is n
and the number of lattice points is m.
If m=n, the conformation of the sequence in the lattice is defined as compact.

In this github repo I define compactness of a connected set of square lattice
points to mean that there are no holes.
